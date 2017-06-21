var panel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style: {width: '300px'}
});

panel.add(ui.Label('Draw features'));

var poly_checkbox = ui.Checkbox({
	label: "Check if project area is available as Google fusion table",
	onChange: function(checked) {
		
	}
});

var poly_ft = ui.Textbox({
	placeholder: 'Enter Google fusion table id',
	onChange: function(text) {
		var poly = ee.String('ft:').cat(text);
	};
});

if(ft_box.getValue()){
	var aoi = ee.FeatureCollection(ft_box.getValue());
}else{
	var aoi = geometry;
}


var slider_cv = ui.Slider({
	min: 0,
	max: 4,
	value: 2.5,
	onChange: function(value){
		var t1 = value;
	}
});

var slider_rcv = ui.Slider({
	min: 0,
	max: 4,
	value: 2.5,
	onChange: function(value){
		var t2 = value;
	}
});

var slider_ndvi = ui.Slider({
	min: 0,
	max: 2,
	value: 0.5,
	onChange: function(value){
		var t3 = value;
	}
});

var slider_nbr = ui.Slider({
	min: 0,
	max: 2,
	value: 0,
	onChange: function(value){
		var t4 = value;
	}
});

var slider_ndsi = ui.Slider({
	min: 0, 
	max: 2, 
	value: 0.5,
	onChange: function(value){
		var t5 = value;
	}
});

panel.add(ft_box);
panel.add(datebox);
panel.add(ui.Label('Use sliders to change default thresholds'));
panel.add(slider_cv);
panel.add(slider_rcv);
panel.add(slider_ndvi);
panel.add(slider_ndsi);

ui.root.add(panel);

//var aoi = range.filterBounds(geometry);
var aoi = range;
//Map.addLayer(change_off.clipToCollection(aoi), {}, 'change');
var today = ee.Date(Date.now());

var bands_re = ee.List(['B1', 'B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'Pan', 'B10', 'TIR1', 'TIR2', 'BQA', 'fmask']);
var bands_orig = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];
var bands_class = ee.List(['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'ndvi', 'ndsi', 'nbr', 'system:time_start']);

var norm = function(img, bnds, region) {
  var imn = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 30,
    maxPixels: 1000175326
  });
  var isd = img.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: region,
    scale:30,
    maxPixels: 1000175326
  });

var array = img.select(bnds).toArray();
imn = imn.toArray(bnds);
isd = isd.toArray(bnds);

return array.subtract(imn).divide(isd).arrayProject([0]).arrayFlatten([bnds]);

};

var ND = function(img, NIR, R, G, SWIR1, SWIR2){
var NBR = img.normalizedDifference([NIR, SWIR2]).rename(["nbr"]);
var NDSI = img.normalizedDifference([G, SWIR1]).rename(["ndsi"]);
//var NNDSI = norm(NDSI, ["ndsi"], region).rename(["nndsi"]);
var NDVI = img.normalizedDifference([NIR, R]).rename(["ndvi"]);
//var NNDVI = norm(NDVI, ["ndvi"], region).rename(["nndvi"]);

return img.addBands(NDSI)
//.addBands(NNDSI)
.addBands(NDVI).addBands(NBR);
//.addBands(NNDVI);
};

var S2maskClouds = function(img) {
  return img.updateMask(img.select('QA60')
  .lt(10)
  //.and(img.select('B9')
  //.lt(1000))
  )
  .addBands(img.metadata('system:time_start')
  );
};

var LSmaskClouds = function(img) {
  var scored = ee.Algorithms.Landsat.simpleCloudScore(img);
  var masked = scored.updateMask(scored.select(['cloud']).lte(20));
  return masked.addBands(img.metadata('system:time_start'));
};

var LSprocess = function(img_col, region, type){
  var area = img_col.filterBounds(region);
  var cloud = area.map(LSmaskClouds);
  var med;
  if (type == 'recent'){
    med = cloud.qualityMosaic('system:time_start')
    .clipToCollection(region);
  } else if (type == 'median'){
      med = cloud.median()
      .clipToCollection(region);
  }
  var rename = med;
  return ND(rename, 'B5', 'B4', 'B3', 'B6', 'B7', aoi);
  //var final = norm(nd, bands_class, aoi);
};

var S2process = function(img_col, region, type){
  var area = img_col.filterBounds(region);
  var cloud = area.map(S2maskClouds);
  var med;
  if (type == 'recent'){
    med = cloud.qualityMosaic('system:time_start')
    .clipToCollection(region);
  } else if (type == 'median'){
      med = cloud.median()
      .clipToCollection(region);
  }
  var rename = med;
  return ND(rename, 'B8', 'B4', 'B3', 'B11', 'B12', aoi);
  //var final = norm(nd, bands_class, aoi);
};

var before_on = S2process(
  S2.filterDate('2015-05-01', '2015-10-01'), aoi, 'median')
  .select(bands_class);
var after_off = S2process(
  S2.filterDate('2016-12-01', '2017-04-01'), aoi, 'median')
  .select(bands_class);
var before_off = S2process(
  S2.filterDate('2015-11-01', '2016-06-01'), aoi, 'median')
  .select(bands_class);
var after_on = S2process(
  S2.filterDate('2016-06-01', '2016-09-30'), aoi, 'median')
  .select(bands_class);

var cv = function(b, a, bnds){
return(
  b.select(bnds)
                  .subtract(a
                            .select(bnds))
                  .pow(2)
.reduce(ee.Reducer.sum())
.rename(['cv'])
);  
};

var rcvmax = function(b, a, bnds){
  return(
    b.select(bnds)
    .subtract(a.select(bnds))
    .divide(b.select(bnds)
     .max(a.select(bnds)))
    .pow(2)
.reduce(ee.Reducer.sum())
.rename(['rcvmax'])
    );
};

var cv_on = cv(before_on, after_on, bands_orig);
var cv_off = cv(before_off, after_off, bands_orig);
var rcv_on = rcvmax(before_on, after_on, bands_orig);
var rcv_off = rcvmax(before_off, after_off, bands_orig);

var d = function(b, a, bnds){
  return(
    b.select(bnds).subtract(a.select(bnds))
    );
};
var diff_on = d(before_on, after_on, ['ndvi', 'nbr', 'ndsi']);
var diff_off = d(before_off, after_off, ['ndvi', 'nbr', 'ndsi']);
var change_on = cv_on.addBands(rcv_on).addBands(diff_on);
var change = cv_off.addBands(rcv_off).addBands(diff_off);

//var test = NLCD.select('landcover').reduceRegions({
//  reducer: ee.Reducer.mode(),
//  collection: non.merge(change),
//  scale: 10,
//});

//var pos_change = image.clipToCollection(geometry).sampleRegions({
//  collection: test,
//  scale: 10
//  });

// Export.table.toDrive({
// collection: pos_change,
// description: 'NRange_chngData',
// fileFormat: 'CSV'
// });

//Export.image.toAsset({image: change_off,
// description: 'Range_S2_LCC_off',
// scale: 10,
// maxPixels:750000000000
//});

var stats = function(img, scl){
  var mean = img.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_mn', 'rcv_mn', 'ndvi_mn', 'ndsi_mn', 'nbr_mn']);

  var sd = img.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_sd', 'rcv_sd', 'ndvi_sd', 'ndsi_sd', 'nbr_sd']);

  var min = img.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: aoi,
  scale: scl,
  maxPixels: 100000000000
  }).rename(['cv', 'rcvmax', 'ndvi', 'ndsi', 'nbr'], ['cv_min', 'rcv_min', 'ndvi_min', 'ndsi_min', 'nbr_min']);  
  
  return mean.combine(sd).combine(min);
};
//print(stats(image, 100));
//print(stats(image2, 100));

var thresh = function(img, t1, t2, t3, t4, t5){
  
  var dic = stats(img, 30);
 var cv_thresh = ee.Number(dic.get('cv_min'))
 .add(ee.Number(dic.get('cv_sd'))
 .multiply(t1));

 var rcv_thresh = ee.Number(dic.get('rcv_min'))
 .add(ee.Number(dic.get('rcv_sd'))
 .multiply(t2));

 var ndvi_thresh = ee.Number(dic.get('ndvi_mn'))
 .add(ee.Number(dic.get('ndvi_sd'))
 .multiply(t3));

 var nbr_thresh = ee.Number(dic.get('nbr_mn'))
 .add(ee.Number(dic.get('nbr_sd'))
 .multiply(t4));
 
var ndsi_thresh = ee.Number(dic.get('ndsi_mn'))
 .add(ee.Number(dic.get('ndsi_sd'))
 .multiply(t5));

 return(
   img.select(['cv']).gt(cv_thresh)
 .and(img.select(['rcvmax']).gt(rcv_thresh))
 .and(img.select(['ndvi']).gt(ndvi_thresh))
 //.and(img.select(['nbr']).gt(nbr_thresh))
 .and(img.select(['ndsi']).lt(ndsi_thresh))
);

};
var thresh2 = function(img, t1, t2, t3, t4, t5){
  
  var dic = stats(img, 30);
 var cv_thresh = ee.Number(dic.get('cv_min'))
 .add(ee.Number(dic.get('cv_sd'))
 .multiply(t1));

 var rcv_thresh = ee.Number(dic.get('rcv_min'))
 .add(ee.Number(dic.get('rcv_sd'))
 .multiply(t2));

 var ndvi_thresh = ee.Number(dic.get('ndvi_mn'))
 .add(ee.Number(dic.get('ndvi_sd'))
 .multiply(t3));

 var nbr_thresh = ee.Number(dic.get('nbr_mn'))
 .add(ee.Number(dic.get('nbr_sd'))
 .multiply(t4));
 
var ndsi_thresh = ee.Number(dic.get('ndsi_mn'))
 .add(ee.Number(dic.get('ndsi_sd'))
 .multiply(t5));

 return(
   img.select(['cv']).gt(cv_thresh)
 .and(img.select(['rcvmax']).gt(rcv_thresh))
 .and(img.select(['ndvi']).gt(ndvi_thresh))
 //.and(img.select(['nbr']).gt(nbr_thresh))
 .and(img.select(['ndsi']).gt(ndsi_thresh))
);
};
var lda = function(img, int, ld1, ld2, ld3, ld4, ldt){
  
  var dic = stats(img, 100);
 var cv_z = img.select('cv').subtract(ee.Number(dic.get('cv_min')))
 .divide(ee.Number(dic.get('cv_sd')))
 .multiply(ld1);

 var rcv_z = img.select('rcvmax').subtract(ee.Number(dic.get('rcv_min')))
 .divide(ee.Number(dic.get('rcv_sd')))
 .multiply(ld2);

 var ndvi_z = img.select('ndvi').subtract(ee.Number(dic.get('ndvi_mn')))
 .divide(ee.Number(dic.get('ndvi_sd')))
 .multiply(ld3);

var ndsi_z = img.select('ndsi').subtract(ee.Number(dic.get('ndvi_mn')))
 .divide(ee.Number(dic.get('ndsi_sd')))
 .multiply(ld4);

var ld = cv_z.add(rcv_z).add(ndvi_z).add(ndsi_z).add(int);

 return(
   ld.gt(ldt)
);

};


var thresh_on = thresh(change_on, t1, t2, t3, 0, t5);
var thresh_off = thresh(change_off.clipToCollection(aoi), 2.5, 2.5, 0.20, 0, -0.5);
//var lda_off = lda(change_off.clipToCollection(geometry), -0.8311436, 0.2370509, -0.2075685, -0.5666503, -0.5299881, 2);
//Map.addLayer(thresh_on, {}, 'thresh_on');
//Map.addLayer(thresh_off, {}, 'thresh_off');
//Map.addLayer(thresh_off.and(thresh_on), {}, 't');
//Map.addLayer(lda_after, {}, 'lda');
//EVALUATION OF THRESHOLDS FOR WELL DETECTION
var validate = function(img){
var cvList = [0.75, 1, 1.5, 2, 2.5];
var rcvList = [0.75, 1, 1.5, 2, 2.5];
var ndviList = [0.25, 0.5, 0.75, 1];
var ndsiList = [0.5, 0.75, 1, 1.5];
var matrix = ee.Dictionary();
for (var i = 0; i < cvList.length; ++i){
 for (var j = 0; j < rcvList.length; ++j){
  for (var k = 0; k < ndviList.length; ++k){
   for (var n = 0; n < ndsiList.length; ++n){
  var index = i.toString().concat(j.toString()).concat(k.toString()).concat(n.toString());
  var thresh_off = thresh(img, cvList[i], rcvList[j], ndviList[k], 0, ndsiList[n]);
  var pos = thresh_off.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: change,
    scale: 10,
    maxPixels:5000000000
  });
  var neg = thresh_off.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: non,
    scale: 10,
    maxPixels:5000000000
  });
  var matrix = matrix.combine(neg.rename(['cv'],['neg'.concat(index)]).combine(pos.rename(['cv'],['pos'.concat(index)])));
   }
  }
 }  
}

return (matrix);
};

var ft = ee.Feature(null, validate(change_off.clipToCollection(geometry)));
var ftCol = ee.FeatureCollection([ft]);

//Export.table.toDrive({
//collection: ftCol,
// description: 'KS_well_evaluation',
// fileFormat: 'CSV'
//});
//Map.addLayer(change, {}, 'change');


//Export.image.toAsset({image: thresh_on,
// description: 'threshold',
// scale: 10,
// maxPixels:750000000000
//});

//Create new attribute with random numbers (0-1]

//Random split for training and testing features
//var training_ft = pads.filterMetadata('split', 'greater_than', 0.25).merge(
// all.filterMetadata('split', 'greater_than', 0.25));

//var training_ft = pads.merge(other);

//var testing_ft = pads.filterMetadata('split', 'less_than', 0.25).merge(
// all.filterMetadata('split', 'less_than', 0.25));

//Sample images for training and testing classifiers

//var training = before.select(bands_class.remove('system:time_start')).sampleRegions(training_ft, ['Class'], 10);

//var testing = all.select(bands).sampleRegions(testing_ft, ['Class'], 10);

//train various classifiers using sample regions and organize into array
//var RF = ee.Classifier.randomForest().train(training, 'Class', bands2); 
//var CNB = ee.Classifier.continuousNaiveBayes().train(training, 'Class', bands2);

//var svm = ee.Classifier.svm().train(training, 'Class', bands_class.remove('system:time_start'));

//var algorithm_list = [RF, CNB, CT];

//define function to return best classifier algorithm
var optimal = function(array){
var classifier;
var trainAccuracy = 0;
//loop to test different classification algorithms, returning the best performer
for (var i = 0; i < array.length; ++i) {
  var validated = training.classify(array[i]);
  var newAccuracy = validated.errorMatrix('Class', 'classification').accuracy();
  if (newAccuracy.gt(trainAccuracy)){
    classifier = array[i];
    trainAccuracy = newAccuracy;
  }
}
return classifier;
};

//return best performing classifier
//var best_classifier = optimal(algorithm_list);
//print(best_classifier);
//classify before and after images
//var class_b = beforeClip.select(bands2).classify(best_classifier);
//var class_a = afterClip.select(bands2).classify(best_classifier);

//var after_svm = after_on.select(bands_class.remove('system:time_start')).classify(svm);
//var before_svm = before.select(bands_class.remove('system:time_start')).classify(svm);

//var after_NB = afterClip.select(bands).classify(CNB);
//var after_RF = afterClip.select(bands).classify(RF);
//var before_CT = beforeClip.select(bands).classify(CT);
//var before_NB = beforeClip.select(bands).classify(CNB);
//var before_RF = beforeClip.select(bands).classify(RF);

//identify changed pixels corresponding to habitat loss
//var change = before_svm.neq(1).and(after_svm.eq(1));
//Map.addLayer(change, {}, 'change');
//Map.addLayer(after_svm, {}, 'after_cl');
//Map.addLayer(before_svm, {}, 'before_cl');

//var opened = lda_off.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'});
//.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'});

var opened2 = thresh_off.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'});
//.focal_mode({radius: 1, kernelType: 'square', units: 'pixels'});

Map.addLayer(opened2, {}, 'thresh');
//Map.addLayer(opened, {}, 'lda');
//var change_poly = opened.eq(1).reduceToVectors({
//  geometry: geometry,
//  scale: 10,
//  geometryType: 'polygon',
//  eightConnected: true,
//  maxPixels: 1000000000
//});
var thresh_poly = opened2.eq(1).reduceToVectors({
  geometry: aoi,
  scale: 10,
  geometryType: 'polygon',
  eightConnected: true,
  maxPixels: 750000000000
});

//Map.addLayer(opened, {}, 'class');

//Map.addLayer(thresh_poly.filter(ee.Filter.eq('label', 1)), {}, 'thresh');
//Export.image.toAsset({
//  image: opened2,
//  description: 'NewWells',
//  scale: 10,
//  maxPixels:750000000000
//});
//Export.table.toDrive(thresh_poly.filter(ee.Filter.eq('label', 1)), 'Well_Confirm_thresh');
//Export.table.toDrive(change_poly.filter(ee.Filter.eq('label', 1)), 'SRange_Wells_lda');