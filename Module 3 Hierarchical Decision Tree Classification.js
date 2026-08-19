//1. Apply Rules from Python Output
// Example rules derived from Python FSMD-CS & SEaTH:
// Layer 1 (Class 1): 7_NDSVI <= 0.1020
// Layer 2 (Class 2): 6_B8A <= 0.2150
// Layer 3 (Class 3): 10_Greenness <= -0.0790
//and so forth.

var unclassified = ee.Image(1); // Mask tracking unclassified pixels (1 = active)
var classMap = ee.Image(0); // Output classification raster initialized to 0

// Layer 1: Extract Class 1 (e.g., Suaeda salsa)
var feat_layer1 = mpsImageMasked.select('7_NDSVI');
var mask_layer1 = feat_layer1.lte(0.1020).and(unclassified);
classMap = classMap.where(mask_layer1, 1);
unclassified = unclassified.and(mask_layer1.not()); // Update unclassified mask

// Layer 2: Extract Class 2 (e.g., Paddy Rice)
var feat_layer2 = mpsImageMasked.select('6_B8A');
var mask_layer2 = feat_layer2.lte(0.2150).and(unclassified);
classMap = classMap.where(mask_layer2, 2);
unclassified = unclassified.and(mask_layer2.not());

// Layer 3: Extract Class 3 (e.g., Maize)
var feat_layer3 = mpsImageMasked.select('10_Greenness');
var mask_layer3 = feat_layer3.lte(-0.0790).and(unclassified);
classMap = classMap.where(mask_layer3, 3);
unclassified = unclassified.and(mask_layer3.not());

//Continue adding layers 4 through N-1 following the same pattern...

// Final Layer: Assign all remaining unclassified pixels to the final class (e.g., Class 8)
classMap = classMap.where(unclassified, 8);

//2. Visualization & Export
var palette = [
  'ff0000', // Class 1
  'ffff00', // Class 2
  '00ff00', // Class 3
  '00ffff', // Class 4
  '0000ff', // Class 5
  'ff00ff', // Class 6
  '800080', // Class 7
  '008000' // Class 8
];

Map.centerObject(roi, 11);
Map.addLayer(classMap.clip(roi), {min: 1, max: 8, palette: palette}, 'Wetland Vegetation Map');

// Export classification map to Google Drive
Export.image.toDrive({
  image: classMap.clip(roi),
  description: 'Wetland_Classification_SEaTH',
  scale: 10,
  region: roi,
  maxPixels: 1e13
});
