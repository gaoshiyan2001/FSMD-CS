//1. Define Parameters & Core Functions
var roi = ee.Geometry.Polygon(/* Your ROI coordinates here */);
var samples = ee.FeatureCollection("users/your_name/wetland_samples"); // Must contain 'class_code' attribute

// Cloud masking function for Sentinel-2 SR
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

// Calculate 25 spectral features (11 bands + 14 spectral indices) for a given month
function addIndices(img, monthStr) {
  var b2 = img.select('B2'), b3 = img.select('B3'), b4 = img.select('B4');
  var b8 = img.select('B8'), b11 = img.select('B11'), b12 = img.select('B12');
  
  var ndvi = img.normalizedDifference(['B8', 'B4']).rename(monthStr + '_NDVI');
  var ndsvi = img.normalizedDifference(['B11', 'B4']).rename(monthStr + '_NDSVI');
  var greenness = img.expression(
    '-0.3599*B2 - 0.3533*B3 - 0.4734*B4 + 0.7243*B8 + 0.0840*B11 - 0.1800*B12',
    {'B2':b2, 'B3':b3, 'B4':b4, 'B8':b8, 'B11':b11, 'B12':b12}
  ).rename(monthStr + '_Greenness');
  
  // Note: Append remaining indices here as needed
  var bands = img.select(['B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12']).regexpRename('^(.*)$', monthStr + '_$1');
                 
  return bands.addBands([ndvi, ndsvi, greenness]);
}

//2. Monthly Composite Generation (March to November)
var months = [3, 4, 5, 6, 7, 8, 9, 10, 11];
var mpsImage = ee.Image([]);

months.forEach(function(m) {
  var startDate = ee.Date.fromYMD(2024, m, 1);
  var endDate = startDate.advance(1, 'month');
  
  var monthlyMedian = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED").filterBounds(roi).filterDate(startDate, endDate).map(maskS2clouds).median();
    
  var monthlyFeatures = addIndices(monthlyMedian, m.toString());
  mpsImage = mpsImage.addBands(monthlyFeatures);
});

// Apply non-vegetation mask (e.g., July NDVI >= 0.2)
var vegMask = mpsImage.select('7_NDVI').gte(0.2);
var mpsImageMasked = mpsImage.updateMask(vegMask);

//3. Sample Region Values & Export CSV
var sampledData = mpsImageMasked.sampleRegions({
  collection: samples,
  properties: ['class_code'],
  scale: 10,
  geometries: false
});

Export.table.toDrive({
  collection: sampledData,
  description: 'MPS_Sample_Features',
  fileFormat: 'CSV'
});
