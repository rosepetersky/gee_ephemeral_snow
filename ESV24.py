##Rose Petersky
##Ephemeral Snow Version 24
##November 8th 2017

#NOTE: Google Earth engine is constantly changing and some functions may be deprecated between 11/8/17 and the day you run this script

import ee
import pandas as pd

# Initialize the Earth Engine object, using the authentication credentials.
ee.Initialize()

##clip imagery to HUC8 watershed boundary
gb = ee.FeatureCollection('ft:1Pdi0p4DK3j9fF7OyFWImLVkZ-DsysG-AAQJMPysO')
gbgeo = gb.geometry()
gabbscoords = gbgeo.coordinates().getInfo()
bigrectangle=gbgeo.bounds()
gbcoords = bigrectangle.coordinates().getInfo()

def gbclip(img):
  gb = ee.FeatureCollection('ft:1Pdi0p4DK3j9fF7OyFWImLVkZ-DsysG-AAQJMPysO')
  gbgeo = gb.geometry()
  imgclip = ee.Image(img).clip(gb)
  return imgclip

##Get the day of the water year
def timeImage(img):
  timestamp = ee.Number(img.get('system:time_start'))
  date = ee.Date(timestamp)
  dayrange = ee.DateRange(ee.Date.fromYMD(year,10,1),ee.Date.fromYMD(wyend,1,1))
  day = ee.Algorithms.If(dayrange.contains(date),ee.Number(date.getRelative('day','year')).double().subtract(ee.Date.fromYMD(year,10,1).difference(ee.Date.fromYMD(year,1,1),'day')),ee.Number(date.getRelative('day','year')).double().add(92))
  doy = ee.Image.constant(day).mask(img.select(0).mask())
  return img.addBands(doy)

##Convert the MODIS bands into a band for snow, a band for no snow, and a band for "other" (clouds,water etc.)
def select_funct(MODimg):
  img = MODimg.select("Fractional_Snow_Cover")
  ##Get all snow and non-snow pixels
  nosnow = img.lt(30)
  snow = img.lte(100).And(img.gte(30))
  invnosnow = nosnow.Not()
  snow = snow.mask(invnosnow)
  ##Region growing method. NOTE: When I ran this, the function I used was ee.Kernel.cross. Now, ee.Kernel.cross is *deprecated* hence the change. 
  snow = snow.reduceNeighborhood(ee.Reducer.allNonZero(),ee.Kernel.rectangle(1,1,"pixels",True,1),"kernel",True,"window")
  ##Get all "other" pixels
  other = img.gt(100)
  snow = snow.unmask(0)
  edges = snow.eq(1).And(other.eq(1))
  ##exclude "region-grown" pixels from the other band
  other = other.where(edges,0)
  ##add all bands together
  bandtest = nosnow.addBands(snow)
  bandtest = bandtest.addBands(other)
  bandtest = bandtest.addBands(MODimg.select("constant"))
  ##bandtest = bandtest.rename(["nosnow","snow","other"])
  bandtest = bandtest.rename(["nosnow","snow","other","doy"])
  return bandtest


##Find all pixels where the snow has appeared or disappeared. Make new bands that show where snow appeared and disappeared each day. 
def find_sndis(img,prev):
  ##prev = ee.List(prev)
  dicty = ee.Dictionary({'collection': ee.List(prev).get(0),'last': ee.List(prev).get(1) })
  ##prev = ee.Dictionary(dicty)
  collection = ee.ImageCollection(dicty.get('collection'))
  ##Get the previous day's image
  last = ee.Image(dicty.get('last'))
  ##Get where snow was present in the previous day's image
  snow = last.select('newsnow')
  ##Get where snow is present in the current day's image
  snow_2 = ee.Image(img).select('snow')
  ##Get where snow is absent in the current day's image 
  nosnow_2 = ee.Image(img).select('nosnow')
  ##Mask all areas where snow is present in the previous day's image
  snow = snow.mask(snow)
  ##Create a band where snow was present in the previous day's image and absent in the current day's image to -1
  sndis = snow.where(nosnow_2,-1)
  ##Change all areas that are equal to -1 to 1 remove the mask on the rest of the areas, and rename the band
  sndis = sndis.eq(-1).unmask(0).rename(['sndis'])
  ##Unmask the "snow" band
  snow = snow.unmask(0)
  ##
  ##Get where snow is absent in the previous day's image 
  nosnow = last.select('newnosnow')
  ##Mask all areas where snow is absent in the previous day's image 
  nosnow = nosnow.mask(nosnow)
  ##Change all areas where snow was absent in the previous day's image and present in the current day's image to -1
  snapp = nosnow.where(snow_2,-1)
  ##Change all areas that are equal to -1 to 1 remove the mask on the rest of the areas, and rename the band
  snapp = snapp.eq(-1).unmask(0).rename(['snapp'])
  ##Unmask the "nosnow" band 
  nosnow = nosnow.unmask(0)
  ##
  ##Get where there are clouds in the current day's image 
  clouds_2 = img.select('other')
  ##Unmask values where there was snow in the previous day's image to 1 and values where there was no snow in the previous day's image to 0
  clouds_2 = ee.Algorithms.If(snow.eq(1),clouds_2.unmask(1),clouds_2.unmask(0))
  ##Create a band where there are clouds in the current day's image AND snow in the previous day's image to 2
  cloudssnow = ee.Image(clouds_2).where(snow,2)
  ##Create a new band where the areas above are equal to 1 and all other areas are equal to 0
  cloudssnow = cloudssnow.eq(2)
  ##Create a band where there are clouds in the current day's image AND no snow in the previous day's image to 2
  cloudsnosnow = ee.Image(clouds_2).where(nosnow,2)
  ##Create a new band where the areas above are equal to 1 and all other areas are equal to 0
  cloudsnosnow = cloudsnosnow.eq(2)
  ##Create a new band that combines the snow in the current image with the areas where there are clouds in the current image and snow in the previous day's image. NOTE: This
  ##method works well for accounting for cloud cover in the Great Basin due to the lack of clouds for most of the region in the spring and summer. Recommend omitting the "cloud snow"
  ##parts of this function if you're in a region with very cloudy springs and summers like the Pacific Northwest and similar. 
  newsnow_2 = snow_2.where(cloudssnow,1).where(nosnow_2,0)
  ##Unmask values where there was snow in the previous day's image to 1 and values where there was no snow in the previous day's image to 0
  newsnow_2 = ee.Algorithms.If(snow.eq(1),newsnow_2.unmask(1),newsnow_2.unmask(0))
  newsnow_2 = ee.Image(newsnow_2).rename(['newsnow'])
  ##
  ##Create a new band that combines the snow in the current image with the areas where there are clouds in the current image and snow in the previous day's image.
  newnosnow_2 = nosnow_2.where(cloudsnosnow,1).where(snow_2,0)
  ##Unmask values where there was no snow in the previous day's image to 1 and values where there was snow in the previous day's image to 0
  newnosnow_2 = ee.Algorithms.If(nosnow.eq(1),newnosnow_2.unmask(1),newnosnow_2.unmask(0))
  newnosnow_2 = ee.Image(newnosnow_2).rename(['newnosnow'])
  ##
  ##Change areas that were flagged for both having snow and no snow into not having snow 
  conflicts = newnosnow_2.eq(1).And(newsnow_2.eq(1))
  newnosnow_2 = newnosnow_2.where(conflicts,0)
  ##
  ##Add bands to the output image collection
  newimg = snapp.addBands(sndis).addBands(newsnow_2).addBands(newnosnow_2)##.addBands(prevsndis).addBands(nosndis)
  newimg = newimg.addBands(ee.Image(img),['nosnow','snow','other','doy'])
  newdicty = ee.Dictionary({'collection': collection.merge(ee.ImageCollection([newimg])),
  'last': newimg})
  listy = newdicty.values(['collection','last'])
  return listy

##Fix errors relating to snow presence at the end of the water year (caps the length of snow covered days at the last day of the water year)
def fixerrors(img,prev):
  dicty = ee.Dictionary({'collection': ee.List(prev).get(0),'last': ee.List(prev).get(1) })
  prev = ee.Dictionary(dicty)
  collection = ee.ImageCollection(dicty.get('collection'))
  last = ee.Image(dicty.get('last'))
  ##
  sndis = ee.Image(last).select('sndis')
  snapp_2 = ee.Image(img).select('snapp')
  prevsndis = sndis.eq(1).And(snapp_2.eq(1)).rename(['prevsndis'])
  ##
  newimg = img.addBands(prevsndis)
  newdicty = ee.Dictionary({'collection': collection.merge(ee.ImageCollection([newimg])),
  'last': newimg})
  listy = newdicty.values(['collection','last'])
  return listy

#Get length of snow covered days for each ephemeral or seasonal event 
def doydiff(img,prev):
  ##prev = ee.List(prev)
  dicty = ee.Dictionary({'collection': ee.List(prev).get(0),'last': ee.List(prev).get(1) })
  ##prev = ee.Dictionary(dicty)
  collection = ee.ImageCollection(dicty.get('collection'))
  last = ee.Image(dicty.get('last'))
  ##select the day of the year
  doy = ee.Image(img).select('doy')
  ##select where snow has disappeared and appeared in this day
  sndis = img.select('sndis').unmask(0)
  ##invert the seleection
  sndisinv = sndis.Not()
  ##change doy values where snow has not disappeared to 0
  doy_sndis = doy.where(sndisinv,0)
  ##select where snow has appeared in this day
  snapp = img.select('snapp').unmask(0)
  ##select everywhere that snow disappeared in the previous image
  prevsndis = img.select('prevsndis')
  ##remove snapp values where snow disappeared in the previous image 
  snapp = snapp.where(prevsndis,0)
  ##invert the selection
  snappinv = snapp.Not()
  ##change doy values where snow has not appeared to 0
  doy_snapp = doy.where(snappinv,0)
  ##select the snow appearance mosaic from the previous day
  prevsnapp = last.select('snappmosaic').unmask(0)
  ##merge the doy values for snow appearance into the previous snow apperance mosaic
  snappmosaic = doy_snapp.add(prevsnapp).unmask(prevsnapp)
  ##subtract the snow appearance mosaic from the snow disappearance doy values
  sndiff = doy_sndis.subtract(snappmosaic)
  ##take all negative difference values (i.e. where snow is still present) and select them
  negs = sndiff.lt(0)
  ##turn the negative values to 0
  sndiff = sndiff.where(negs,0).rename(['sndiff'])
  ##saves the current snappmosaic
  snappmosaic_2 = snappmosaic.rename(['snappmosaic_2'])
  prevsnappmosaic_2 = last.select('snappmosaic_2')
  ##saves the snappmosaic instead of resetting values where snow disappeared to zero (to get day of first snow appearance) 
  firstsnow = snappmosaic
  ##set the snow appearence doy values to 0 everywhere that snow has disappeared 
  snappmosaic = snappmosaic.where(sndis,0).rename(['snappmosaic'])
  ##changes the snappmosaic everywhere snow disappeared in the previous image to snappmosaic_2
  snappmosaic = snappmosaic.where(prevsndis,prevsnappmosaic_2)
  newimg = img.addBands(snappmosaic).addBands(snappmosaic_2)
  ##invert the prevsndis band
  prevsndisinv = prevsndis.Not()
  ##get the previous sndiff values
  prevsndiff = last.select('sndiff')
  ##use the inverted band to get the previous sndiff values for all erroroneous events
  badsndiff = prevsndiff.where(prevsndisinv,0).rename(['badsndiff'])
  ##add the bad sndiff to the image
  newimg=newimg.addBands(badsndiff)
  ##add the new sndiff to the image
  newimg = newimg.addBands(sndiff)
  ##add firstsnow to the image
  newimg.addBands(firstsnow).rename(["snapp","sndis","newsnow","newnosnow","nosnow", "snow", "other", "doy","prevsndis","snappmosaic","snappmosaic_2","badsndiff","sndiff","firstsnow"])
  newdicty = ee.Dictionary({'collection': collection.merge(ee.ImageCollection([newimg])),
  'last': newimg})
  listy = newdicty.values(['collection','last'])
  return listy

#Remove all ephemeral events with a length of only one day (removed due to prevalence of anomalies that affect one day ephemeral events) 
def goodsndiffs(img,prev):
  ##prev = ee.List(prev)
  dicty = ee.Dictionary({'collection': ee.List(prev).get(0),'last': ee.List(prev).get(1) })
  ##prev = ee.Dictionary(dicty)
  collection = ee.ImageCollection(dicty.get('collection'))
  last = ee.Image(dicty.get('last'))
  prevsndiff=last.select('sndiff')
  badsndiff=ee.Image(img).select('badsndiff')
  goodsndiff=prevsndiff.where(badsndiff,0).rename(['goodsndiff']).unmask(0)
  #
  newimg=img.addBands(goodsndiff)
  newdicty = ee.Dictionary({'collection': collection.merge(ee.ImageCollection([newimg])),
  'last': newimg})
  listy = newdicty.values(['collection','last'])
  return listy

#Count the number of snow covered days (used to get average SSM) 
def sndiff_count(img,prev):
  dicty = ee.Dictionary({'collection': ee.List(prev).get(0),'last': ee.List(prev).get(1) })
  ##prev = ee.Dictionary(dicty)
  collection = ee.ImageCollection(dicty.get('collection'))
  last = ee.Image(dicty.get('last'))
  #
  goodsndiff=ee.Image(img).select('goodsndiff')#.unmask(0)
  count=goodsndiff.gt(0)
  prevcount=last.select('countband')
  countband=count.add(prevcount).rename(['countband'])
  #
  newimg=img.addBands(countband)
  newdicty = ee.Dictionary({'collection': collection.merge(ee.ImageCollection([newimg])),
  'last': newimg})
  listy = newdicty.values(['collection','last'])
  return listy

#Get the number of seasonal and ephemeral snow events
def event_count(img):
  goodsndiff = ee.Image(img).select('goodsndiff')
  #badsndiff = ee.Image(img).select('badsndiff')
  e_events = goodsndiff.gte(1).And(goodsndiff.lt(60))
  s_events = goodsndiff.gte(60)
  zeros = ee.Image.constant(0)
  #Get the SSM numerator value (this value is positive during a seasonal snow event and negative during an ephemeral snow event)
  snowscale = goodsndiff.where(e_events,zeros.subtract(goodsndiff))
  #bad_e_events = badsndiff.gte(1).And(badsndiff.lte(60))
  #bad_s_events = badsndiff.gt(60)
  newimg = e_events.addBands(s_events).addBands(snowscale)
  newimg = newimg.rename(['e_events','s_events','snowscale'])
  return newimg

#SNOTEL stations in the Great Basin
sntlstations=ee.FeatureCollection('ft:1cpG0rQqZYLMAcawEY4SBbYrZIMFxHK2jpi8R1N1U')

#This can be changed to any year between 2000 and the present water year
yearlist = [2015]

#stnlist=pd.read_csv('/Users/rpetersky/Documents/ES_Project/stnlist.csv')

#print(stnlist)

for year in yearlist:
  #Get the value matching the end of the water year
  wyend=year+1

  #Get the MODIS images for the water year
  MODIS = ee.ImageCollection("MODIS/MOD10A1").filterDate(str(year)+'-10-01', str(wyend)+'-09-30')

    ##list = ee.ImageCollection(MODIS).toList(300)

    ##img2 = list.get(35)

  #Clip the MODIS images to the Great Basin
  MODclip = MODIS.map(gbclip)

  #Get the day of the water year for each image
  MODISdoy = MODclip.map(timeImage)

  #Get the snow, no snow, and other bands 
  multicol = MODISdoy.map(select_funct)

    ##list = ee.ImageCollection(multicol).toList(175)

    ##img2 = list.get(84)

  #Get snow appearence and disappearence 
  multicol_sd = multicol.iterate(find_sndis, ee.Dictionary({'collection': ee.ImageCollection([]),
    'last': ee.Image([0, 0, 0, 0, 0, 0, 0, 0]).select([".*"], 
    ["snapp","sndis","newsnow","newnosnow","nosnow", "snow", "other", "doy"])  }).values(['collection','last']))
  
  sd = ee.List(multicol_sd).get(0)

    ##list = ee.ImageCollection(sd).toList(300)

    ##img3 = list.get(84)

    ##img3_nosnow = ee.Image(img3).select['newsnow']

  #Fix errors related to snow appearance and disappearance at the end of the water year
  getprevsndis = ee.ImageCollection(sd).iterate(fixerrors, ee.Dictionary({'collection': ee.ImageCollection([]),
    'last': ee.Image([0, 0, 0, 0, 0, 0, 0, 0, 0]).select([".*"], 
    ["snapp","sndis","newsnow","newnosnow","nosnow", "snow", "other", "doy","prevsndis"])  }).values(['collection','last']))

  withprevsndis = ee.List(getprevsndis).get(0)

  #Get the lengths of snow covered days 
  sndis_events = ee.ImageCollection(withprevsndis).iterate(doydiff, ee.Dictionary({'collection': ee.ImageCollection([]),
    'last': ee.Image([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).select([".*"], 
    ["snapp","sndis","newsnow","newnosnow","nosnow", "snow", "other", "doy","prevsndis","snappmosaic","snappmosaic_2","badsndiff","sndiff","firstsnow"])  }).values(['collection','last']))

  sndiffs_col = ee.List(sndis_events).get(0)

  #Remove all ephemeral events that are 1 day long
  getgoodsndiff = ee.ImageCollection(sndiffs_col).iterate(goodsndiffs, ee.Dictionary({'collection': ee.ImageCollection([]),
    'last': ee.Image([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]).select([".*"], 
    ["snapp","sndis","newsnow","newnosnow","nosnow", "snow", "other", "doy","prevsndis","snappmosaic","snappmosaic_2","badsndiff","sndiff","firstsnow","goodsndiff"])  }).values(['collection','last']))

  goodsndiffcol=ee.List(getgoodsndiff).get(0)

  goodsndiff = ee.ImageCollection(goodsndiffcol).select('goodsndiff')

  #Get the maximum number of snow covered days 
  maxsndiff = goodsndiff.cast({'goodsndiff':'float'},['goodsndiff']).max()

  #Get the date of first snow 
  firstsnow = ee.ImageCollection(goodsndiffcol).select('firstsnow').cast({'firstsnow':'float'},['firstsnow']).max()

  #Get the SSM, number of ephemeral snow events, and number of seasonal snow events 
  ephemeral_func = ee.ImageCollection(goodsndiffcol).map(event_count)

  #Create images contaning the total number of ephemeral and seasonal snow events 
  ephemeral = ee.ImageCollection(ephemeral_func).select('e_events').sum().cast({'e_events':'float'})
  seasonal = ee.ImageCollection(ephemeral_func).select('s_events').sum().cast({'s_events':'float'})

  #Get the total number of snow covered days
  sndiffsum = ee.ImageCollection(sndiffcol).select('goodsndiff').cast({'goodsndiff':'float'},['goodsndiff']).sum()

  #Create an image for the SSM (numerator is the 'snowscale' band and the denominator is the total number of snow covered days) 
  snowscale = ee.ImageCollection(ephemeral_func).select('snowscale').cast({'snowscale':'float'},['snowscale']).sum().divide(sndiffsum)

  #Export maximum snow covered days as an image. This can be changed to the total number of ephemeral events, seasonal events, or the SSM by changing the image that goes in the first
  #spot in ee.batch.Export.image(). Change the name by changing the second entry. Change the scale (pixels), maximum pixels, or region by changing the respective entries in the dictionary. 
  task = ee.batch.Export.image(maxsndiff,'maxsndiff_'+str(wyend)+'_11_30',{'scale':500,'maxPixels':1e12,'region':str(gbcoords)})
  task.start() 
  while task.status()['state'] == 'READY':
    print 'Running...'




















