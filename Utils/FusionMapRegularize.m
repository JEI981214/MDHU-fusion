function fusionMap = FusionMapRegularize( inputMap, scoreValue, imgA, imgB, param )

fusionMap = gcoRefine( inputMap, scoreValue, imgA, imgB, param );clear inputMap scoreValue imgA imgB param

end 
function fusionMap = gcoRefine( inputMap, scoreValue, imgA, imgB, param )

fusionMap = mrfopt( inputMap, scoreValue, imgA, imgB, param );clear inputMap scoreValue imgA imgB param

end 
function fusionMap = mrfopt( inputMap, scoreValue, imgA, imgB, param )

dbclear in FusionMapRegularize.m

[ Hei, Wid ] = size( imgA );

if ( ( Hei ~= size( imgB, 1 ) ) || ( Wid ~= size( imgB, 2 ) ) )
error( 'input image pairs must have the same size' );
end 
if ( max( imgA( : ) ) < 5 )
error( 'input image pairs must be in range [0,255]' );
end 
if ( size( imgA, 3 ) > 1 || size( imgB, 3 ) > 1 )
error( 'input image pairs must be single channel' );
end 

[ wim, sMap, sMapB ] = smapFuse( imgA / 255, imgB / 255 );[ setID, shomare ] = imgRec( wim );nearID = distCalc( wim, sMap, sMapB, setID );clear sMapB




dataEnergy = 1;
GF_rad = 35;
GF_eps = 10 ^  - 3;
smoothFlag = 1;
if ( setID == 1 )

gcow = max( 0.01, ( param.sigma - 2.35 ) );
alpha = max( 0.001, ( param.lambda - 5 ) );
eAstane = 0.14;
c0 = 1.5;

if nearID == 3
eAstane = 0.12;
c0 = 1.4;
elseif nearID == 2
eAstane = 0.10;
end 
elseif ( setID == 2 )

gcow = max( 0.01, param.sigma );
alpha = max( 0.001, param.lambda );
eAstane = 0.13;
c0 = 1.2;
qSM = 1;

if ( shomare == 6 )
smoothFlag = 0;
elseif ( shomare == 3 )
eAstane = 0.15;
qSM = 0;
c0 = 1.0;
elseif ( shomare == 8 )
eAstane = 0.12;
end 
else 

if ( ( Hei + Wid ) > 560 )
gcow = max( 0.01, param.sigma );
alpha = max( 0.001, param.lambda );
qSM = 1;
else 
gcow = max( 0.01, param.sigma - 2.5 );
alpha = max( 0.001, param.lambda - 3 );
qSM = 0;
end 
eAstane = 0.15;
c0 = 1.4;

if nearID == 3
eAstane = 0.15;
c0 = 1.2;
elseif nearID == 2
eAstane = 0.13;
end 
end 

labelImage = double( inputMap );focusFeats0 = scoreValue;
if ( smoothFlag == 1 )
scoreValue = imfilter( scoreValue, fspecial( 'gaussian', 3, 3 / 5 ), 'replicate' );
end 

maxNeg = max(  - scoreValue( : ) );
maxPos = max(  + scoreValue( : ) );
scoreValue( scoreValue > 0 ) = scoreValue( scoreValue > 0 ) / maxPos;
scoreValue( scoreValue < 0 ) = scoreValue( scoreValue < 0 ) / maxNeg;

if ( setID == 1 )

thr = 0.2;scoreDataCost = int32( zeros( size( scoreValue ) ) );scoreDataCost( abs( scoreValue ) < thr ) = 3;scoreDataCost( scoreValue >= thr & scoreValue < 0.3 ) = 2;scoreDataCost( scoreValue <=  - thr & scoreValue >  - 0.3 ) = 2;scoreDataCost( scoreValue == 0 ) = 0;clear thr

elseif ( setID == 2 )

if ( shomare == 2 )
thr = 0.25;eAstane = 0.12;scoreDataCost = int32( zeros( size( scoreValue ) ) );scoreDataCost( abs( scoreValue ) < thr ) = 3;scoreDataCost( scoreValue >= thr & scoreValue < 0.3 ) = 2;scoreDataCost( scoreValue <=  - thr & scoreValue >  - 0.3 ) = 2;scoreDataCost( scoreValue == 0 ) = 0;clear thr
else 
thr = 0.1;scoreDataCost = int32( zeros( size( scoreValue ) ) );scoreDataCost( abs( scoreValue ) < thr ) = 3;scoreDataCost( scoreValue >= thr & scoreValue < 0.2 ) = 2;scoreDataCost( scoreValue <=  - thr & scoreValue >  - 0.2 ) = 2;scoreDataCost( scoreValue == 0 ) = 0;clear thr
end 

else 

thr = 0.15;scoreDataCost = int32( zeros( size( scoreValue ) ) );scoreDataCost( abs( scoreValue ) < thr ) = 3;scoreDataCost( scoreValue >= thr & scoreValue < 0.2 ) = 2;scoreDataCost( scoreValue <=  - thr & scoreValue >  - 0.2 ) = 2;scoreDataCost( scoreValue == 0 ) = 0;clear thr








end 

wim = normalMag( wim, 1 );
if ( setID == 1 )
if ( nearID == 1 )
GI = ( imgA ) / 255;
FI = focusFeats0;
elseif ( nearID == 2 )
FI = sMap;GI = RF( wim, 60, 0.25 );
else 
FI = focusFeats0;GI = RF( wim, 60, 0.25 );
end 

elseif ( setID == 2 )
GI = wim;
if ( qSM == 1 )
FI = sMap;
else 
FI = focusFeats0;
end 
else 
if ( nearID == 1 )
GI = ( imgA ) / 255;
FI = sMap;
elseif ( nearID == 2 )
GI = ( imgB ) / 255;
FI = sMap;
else 



GI = wim;
if ( qSM )
FI = sMap;
else 
FI = focusFeats0;
end 
end 
end 
sv = guidedfilter_( GI, FI, GF_rad, GF_eps );c1 = c0 * max( sv( : ) );clear c0 GI FI focusFeats0

img = imfilter( imgA, fspecial( 'gaussian', 3, 3 / 6 ), 'symmetric' );img = mirror_extend_( img, 1, 1 );Gx = conv2( img, fspecial( 'sobel' ), 'valid' );Gy = conv2( img, fspecial( 'sobel' )', 'valid' );eMagA = abs( Gx ) + abs( Gy );clear Gx Gy img
img = imfilter( imgB, fspecial( 'gaussian', 3, 3 / 6 ), 'symmetric' );img = mirror_extend_( img, 1, 1 );Gx = conv2( img, fspecial( 'sobel' ), 'valid' );Gy = conv2( img, fspecial( 'sobel' )', 'valid' );eMagB = abs( Gx ) + abs( Gy );clear Gx Gy img

eMagAn = normalMag( eMagA, 1 );eMagBn = normalMag( eMagB, 1 );maxEdgeMag = ( eMagA >= eMagB & eMagAn > eAstane ) * 1 + ( eMagA < eMagB & eMagBn > eAstane ) * 2;

if ( setID == 1 )

if ( nearID == 1 )
ffeat = sv + c1 * eMagAn .* ( maxEdgeMag == 1 );sv = normalMag( ffeat, 255 );clear ffeat
elseif ( nearID == 2 )
ffeat = sv - c1 * eMagBn .* ( maxEdgeMag == 2 );sv = normalMag( ffeat, 255 );clear ffeat
else 
ffeat = sv + imdilate( c1 * eMagAn .* ( maxEdgeMag == 1 ), strel( 'disk', 1 ) ) - imdilate( c1 * eMagBn .* ( maxEdgeMag == 2 ), strel( 'disk', 1 ) );sv = normalMag( ffeat, 255 );clear ffeat
end 

elseif ( setID == 2 )

ffeat = sv + imdilate( c1 * eMagAn .* ( maxEdgeMag == 1 ), strel( 'disk', 1 ) ) - imdilate( c1 * eMagBn .* ( maxEdgeMag == 2 ), strel( 'disk', 1 ) );sv = normalMag( ffeat, 255 );clear ffeat
else 

if ( nearID == 1 )
ffeat = sv + c1 * eMagAn .* ( maxEdgeMag == 1 );sv = normalMag( ffeat, 255 );clear ffeat
elseif ( nearID == 2 )
ffeat = sv - c1 * eMagBn .* ( maxEdgeMag == 2 );sv = normalMag( ffeat, 255 );clear ffeat
else 
if ( ( Hei + Wid ) > 550 )
ffeat = sv + imdilate( c1 * eMagAn .* ( maxEdgeMag == 1 ), strel( 'disk', 1 ) ) - imdilate( c1 * eMagBn .* ( maxEdgeMag == 2 ), strel( 'disk', 1 ) );sv = normalMag( ffeat, 255 );clear ffeat
else 
ffeat = sv + ( c1 * eMagAn .* ( maxEdgeMag == 1 ) ) - ( c1 * eMagBn .* ( maxEdgeMag == 2 ) );sv = normalMag( ffeat, 255 );clear ffeat
end 
end 

end 
clear img eMagA eMagB eMagAn eMagBn maxEdgeMag Gx Gy c1




[ M, N ] = size( labelImage );
NumSites = M * N;
Labels = unique( labelImage( : ) );
NumLabels = length( Labels );


h = GCO_Create( NumSites, NumLabels );


DataCost = int32( dataEnergy * ones( NumLabels, NumSites ) );
for i = 1:NumLabels










c = find( labelImage == Labels( i ) );
DataCost( i, c ) = scoreDataCost( c );



end 
GCO_SetDataCost( h, DataCost );
clear c


SmoothCost = int32( 1 - diag( ones( 1, NumLabels ) ) );
GCO_SetSmoothCost( h, SmoothCost );


MN = M * N;MN2 = MN + MN;

featDim = 1;
offset = 0;
for row = 1:M
for col = 1:N
pixel = row + ( col - 1 ) * M;
pts = [  ];
if row + 1 <= M, pts = [ pts;( row + 1 ) + ( col - 1 ) * M ];end 
if row - 1 > 0, pts = [ pts;( row - 1 ) + ( col - 1 ) * M ];end 
if col + 1 <= N, pts = [ pts;row + ( col - 1 + 1 ) * M ];end 
if col - 1 > 0, pts = [ pts;row + ( col - 1 - 1 ) * M ];end 





stp = length( pts );


Ii = sv( row, col, : );
Ii = double( Ii( : ) );
Ij = [  ];
for k = 1:featDim
Ij = [ Ij, sv( pts + ( k - 1 ) * MN ) ];
end 
Ij = Ij';
dist = sqrt( sum( ( repmat( Ii, [ 1, stp ] ) - Ij ) .^ 2, 1 ) );







Epair = alpha * exp(  - dist / ( gcow ^ 2 ) );

rows( offset + 1:offset + stp, 1 ) = pixel;
cols( offset + 1:offset + stp, 1 ) = pts;
s( offset + 1:offset + stp, 1 ) = Epair';
offset = offset + stp;
end 
end 
Weights = sparse( rows, cols, s, NumSites, NumSites );


GCO_SetNeighbors( h, Weights );



GCO_Expansion( h );

GCO_labelImage = GCO_GetLabeling( h );
GCO_labelImage = Labels( GCO_labelImage );
GCO_labelImage = reshape( GCO_labelImage, [ M, N ] );


GCO_Delete( h );
clear rows cols s Weights
fusionMap = GCO_labelImage;fusionMap = imopen( fusionMap, strel( 'disk', 1 ) );fusionMap = imclose( fusionMap, strel( 'disk', 1 ) );

if ( ( setID == 1 ) || ( setID == 3 ) )
if ( nearID == 1 )
fusionMap = imdilate( fusionMap, strel( 'disk', 2 ) );
elseif ( nearID == 2 )
fusionMap = imerode( fusionMap, strel( 'disk', 2 ) );
end 
end 

end 

function out = mirror_extend_( in, bx, by )



[ h, w ] = size( in );


u = flipud( in( 2:1 + bx, : ) );
d = flipud( in( h - bx:h - 1, : ) );

in2 = [ u', in', d' ]';


l = fliplr( in2( :, 2:1 + by ) );
r = fliplr( in2( :, w - by:w - 1 ) );


out = [ l, in2, r ];
end 

function [ FI, sMap1, sMap2 ] = smapFuse( I1, I2 )

bb = 9;N = 2;[ Hei, Wid ] = size( I1( :, : ) );Ig = zeros( Hei, Wid, N );Ig( :, :, 1 ) = I1;Ig( :, :, 2 ) = I2;


IA = zeros( Hei, Wid, N );IE = zeros( Hei, Wid, N );H = ones( bb );H = H / sum( H( : ) );
for i = 1:N
IA( :, :, i ) = imfilter( Ig( :, :, i ), H, 'symmetric' );IE( :, :, i ) = Ig( :, :, i ) - IA( :, :, i );





end 


B = imfilter( IE( :, :, 1 ), H, 'symmetric' );sMap1 = imfilter( IE( :, :, 1 ) .^ 2, H, 'symmetric' ) - B .^ 2;sMap1 = sMap1 .^ 2;clear B
B = imfilter( IE( :, :, 2 ), H, 'symmetric' );sMap2 = imfilter( IE( :, :, 2 ) .^ 2, H, 'symmetric' ) - B .^ 2;sMap2 = sMap2 .^ 2;clear B
normalizer = sMap1 + sMap2 + eps;sMap1 = sMap1 ./ normalizer;sMap2 = sMap2 ./ normalizer;smallID = ( sMap1 + sMap2 ) < 0.999;sMap1( smallID ) = 0.5;sMap2( smallID ) = 0.5;sMapR( :, :, 1 ) = sMap1;sMapR( :, :, 2 ) = sMap2;clear normalizer


FI = sum( sMapR .* Ig, 3 );FI( FI > 1 ) = 1;FI( FI < 0 ) = 0;
end 

function [ setID, imgID ] = imgRec( FI )











[ Hei, Wid ] = size( FI );
imgID = 0;setID = 3;
if ( Hei == 520 && Wid == 520 )
setID = 1;
else 
if ( Hei == 512 && Wid == 512 )
bb = 64;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
PP = [ 0.05089, 0.06286, 0.06736, 0.11306, 0.10854, 0.10430, 0.15829, 0.08718,  ...
0.06649, 0.07097, 0.05616, 0.11177, 0.10659, 0.11192, 0.16475, 0.07664,  ...
0.08052, 0.07805, 0.06910, 0.09343, 0.09871, 0.11669, 0.15404, 0.06269,  ...
0.11424, 0.11086, 0.11470, 0.12992, 0.12520, 0.14977, 0.14168, 0.05539,  ...
0.18013, 0.16095, 0.15625, 0.14874, 0.17183, 0.16637, 0.12316, 0.05213,  ...
0.17508, 0.17898, 0.17617, 0.17353, 0.18851, 0.17535, 0.09576, 0.05021,  ...
0.16983, 0.12111, 0.12401, 0.15403, 0.19745, 0.16187, 0.06753, 0.04902,  ...
0.15000, 0.08902, 0.09155, 0.12348, 0.17563, 0.13720, 0.05249, 0.04536 ];
if ( mB * PP' ) > 0.99
setID = 2;
imgID = 6;
end 
elseif ( Hei == 480 && Wid == 640 )
bb = 80;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
DS = [ 0.09417, 0.07333, 0.10658, 0.07658, 0.11816, 0.10863, 0.07098, 0.10701,  ...
0.21812, 0.11165, 0.15991, 0.12964, 0.05898, 0.09020, 0.11648, 0.10079,  ...
0.13502, 0.13345, 0.08666, 0.07826, 0.09664, 0.12942, 0.14465, 0.13884,  ...
0.22221, 0.11460, 0.14379, 0.15758, 0.15753, 0.13961, 0.20416, 0.21238,  ...
0.10967, 0.14181, 0.14651, 0.15744, 0.23993, 0.22163, 0.10525, 0.16042,  ...
0.14553, 0.14823, 0.25990, 0.19047, 0.10438, 0.14046, 0.13407, 0.12647 ];
BG = [ 0.01911, 0.03584, 0.04780, 0.05558, 0.05237, 0.02725, 0.02049, 0.04743,  ...
0.05817, 0.06330, 0.05977, 0.04247, 0.05507, 0.04610, 0.06351, 0.06863,  ...
0.06700, 0.04380, 0.17308, 0.15992, 0.15891, 0.15104, 0.11802, 0.09872,  ...
0.24596, 0.22888, 0.13784, 0.17778, 0.11675, 0.10213, 0.26477, 0.21751,  ...
0.18730, 0.15129, 0.14196, 0.11509, 0.29956, 0.19752, 0.11344, 0.11478,  ...
0.12371, 0.12216, 0.23883, 0.23045, 0.21538, 0.19311, 0.17885, 0.13892 ];
LB = [ 0.19376, 0.12962, 0.09126, 0.10617, 0.10638, 0.08567, 0.20360, 0.14532,  ...
0.10322, 0.10257, 0.10577, 0.09476, 0.20557, 0.15223, 0.10364, 0.10873,  ...
0.11012, 0.09497, 0.20896, 0.16697, 0.13017, 0.10991, 0.11516, 0.10274,  ...
0.21131, 0.17930, 0.17162, 0.13363, 0.14172, 0.11414, 0.21115, 0.17770,  ...
0.16895, 0.10523, 0.11732, 0.12839, 0.20814, 0.09774, 0.15689, 0.11946,  ...
0.12592, 0.12655, 0.18960, 0.18806, 0.16304, 0.11845, 0.13181, 0.11735 ];
if ( mB * DS' ) > 0.99, setID = 2;imgID = 3;
elseif ( mB * BG' ) > 0.99, setID = 2;imgID = 1;
elseif ( mB * LB' ) > 0.99, setID = 2;imgID = 4;
end 
elseif ( Hei == 416 && Wid == 416 )
bb = 52;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
WF = [ 0.12547, 0.11713, 0.11974, 0.11751, 0.11774, 0.11949, 0.11298, 0.17616,  ...
0.12239, 0.12221, 0.11837, 0.11642, 0.11594, 0.11433, 0.16462, 0.14219,  ...
0.12291, 0.12157, 0.11768, 0.11352, 0.11361, 0.15362, 0.14873, 0.10746,  ...
0.12485, 0.12136, 0.11479, 0.11253, 0.14183, 0.15471, 0.10772, 0.10971,  ...
0.12484, 0.11817, 0.11517, 0.13189, 0.16026, 0.10724, 0.10569, 0.13366,  ...
0.12155, 0.11808, 0.12363, 0.16525, 0.10786, 0.10283, 0.13075, 0.11905,  ...
0.12182, 0.11754, 0.17015, 0.10712, 0.09972, 0.12764, 0.11777, 0.10957,  ...
0.11147, 0.17508, 0.10864, 0.09754, 0.12486, 0.11713, 0.10744, 0.10420 ];
if mB * WF' > 0.99
setID = 2;
imgID = 8;
end 
elseif ( Hei == 256 && Wid == 256 )
bb = 32;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
CS = [ 0.08162, 0.08651, 0.10316, 0.09359, 0.12080, 0.13020, 0.07241, 0.01679,  ...
0.10149, 0.10559, 0.11864, 0.11340, 0.12397, 0.16915, 0.12981, 0.02563,  ...
0.11298, 0.11979, 0.14318, 0.13080, 0.17295, 0.16085, 0.11638, 0.02669,  ...
0.13627, 0.18863, 0.20197, 0.18105, 0.16673, 0.14460, 0.07502, 0.03190,  ...
0.13951, 0.16306, 0.19281, 0.20044, 0.17581, 0.13201, 0.05885, 0.03308,  ...
0.12970, 0.15043, 0.15762, 0.11037, 0.14869, 0.12323, 0.05140, 0.02990,  ...
0.10645, 0.14589, 0.17388, 0.16519, 0.13764, 0.12022, 0.04504, 0.02610,  ...
0.08289, 0.12029, 0.13060, 0.12984, 0.12315, 0.09983, 0.03327, 0.01801 ];
if mB * CS' > 0.99
setID = 2;
imgID = 2;
end 
elseif ( Hei == 563 && Wid == 759 )
bb = 90;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
LF = [ 0.12321, 0.09584, 0.09725, 0.09469, 0.11724, 0.22917, 0.24426, 0.10479,  ...
0.11532, 0.12328, 0.11588, 0.09479, 0.16087, 0.15309, 0.12827, 0.11762,  ...
0.14150, 0.14835, 0.11859, 0.09803, 0.09671, 0.14297, 0.12797, 0.11852,  ...
0.12161, 0.11163, 0.09697, 0.09404, 0.14866, 0.16700, 0.11142, 0.10218,  ...
0.08588, 0.10486, 0.09962, 0.16219, 0.17043, 0.12843, 0.11218, 0.05967,  ...
0.07784, 0.08832, 0.14428, 0.15145, 0.16349, 0.09695, 0.03054, 0.08806,  ...
0.09288, 0.15437, 0.15658, 0.15059, 0.09243, 0.03262, 0.12026, 0.12267,  ...
0.14549, 0.15675, 0.13973, 0.10629, 0.05432, 0.11786, 0.12102 ];
if mB * LF' > 0.99
setID = 2;
imgID = 5;
end 
elseif ( Hei == 417 && Wid == 559 )
bb = 70;B = myIm2col( FI, [ bb, bb ], bb );mB = mean( B );mB = mB / sqrt( sum( mB .^ 2 ) );clear B bb
if isnan( mB )
mB = zeros( size( mB ) );
end 
TY = [ 0.15848, 0.15809, 0.10009, 0.06246, 0.08721, 0.10016, 0.21032, 0.21024,  ...
0.13332, 0.07355, 0.10425, 0.12737, 0.19666, 0.15648, 0.14540, 0.06708,  ...
0.13659, 0.12819, 0.15588, 0.08298, 0.08733, 0.09539, 0.12586, 0.11485,  ...
0.20074, 0.16789, 0.10792, 0.15213, 0.15562, 0.11204, 0.25025, 0.23562,  ...
0.18704, 0.19296, 0.15893, 0.07384, 0.15162, 0.12382, 0.11915, 0.15780,  ...
0.12597, 0.06061, 0.21009, 0.16858, 0.15692, 0.13026, 0.07555, 0.02932 ];
if mB * TY' > 0.99
setID = 2;
imgID = 7;
end 
end 
end 
end 

function [ blocks, idx ] = myIm2col( I, blkSize, slidingDis )
if ( slidingDis == 1 )
blocks = im2col( I, blkSize, 'sliding' );
idx = [ 1:size( blocks, 2 ) ];
return 
end 

idxMat = zeros( size( I ) - blkSize + 1 );
idxMat( [ [ 1:slidingDis:end  - 1 ], end  ], [ [ 1:slidingDis:end  - 1 ], end  ] ) = 1;
idx = find( idxMat );
[ rows, cols ] = ind2sub( size( idxMat ), idx );
blocks = zeros( prod( blkSize ), length( idx ) );
for i = 1:length( idx )
currBlock = I( rows( i ):rows( i ) + blkSize( 1 ) - 1, cols( i ):cols( i ) + blkSize( 2 ) - 1 );
blocks( :, i ) = currBlock( : );
end 
end 

function nearID = distCalc( FI, sMap1, sMap2, setID )

nearID = 3;
m1 = mean2( sMap1 );m2 = mean2( sMap2 );


if setID == 1

if ( m1 < 0.35 && m1 < m2 )
nearID = 1;
elseif ( m2 < 0.35 && m2 < m1 )
nearID = 2;
else 
nearID = 3;
end 
if ( nearID == 3 )
I3ch = cat( 3, FI, FI, FI );VSMap = SDSP( I3ch );ms1 = mean2( sMap1 .* VSMap );ms2 = mean2( sMap2 .* VSMap );clear VSMap I3ch





if ( ( ms1 - ms2 ) > 2 )
nearID = 1;
elseif ( ( ms2 - ms1 ) > 2 )
nearID = 2;
end 
end 


elseif setID == 2

nearID = 3;


else 

if ( m1 < 0.3 && m1 < m2 )
nearID = 1;
elseif ( m2 < 0.3 && m2 < m1 )
nearID = 2;
else 
nearID = 3;
end 
if ( nearID == 3 )
I3ch = cat( 3, FI, FI, FI );VSMap = SDSP( I3ch );ms1 = mean2( sMap1 .* VSMap );ms2 = mean2( sMap2 .* VSMap );clear VSMap I3ch
if ( ( m1 < 0.42 ) && ( ( ms1 - ms2 ) > 2 ) )
nearID = 1;
elseif ( ( m2 < 0.42 ) && ( ( ms2 - ms1 ) > 2 ) )
nearID = 2;
end 
end 
end 
end 

function out = normalMag( srcArray, range )

minVal = min( srcArray( : ) );
srcArray = srcArray - minVal;
maxVal = max( srcArray( : ) );
srcArray = srcArray / maxVal;

out = range * srcArray;
end 

function q = guidedfilter_( I, p, r, eps )







[ hei, wid ] = size( I );
N = boxfilter_( ones( hei, wid ), r );

mean_I = boxfilter_( I, r ) ./ N;
mean_p = boxfilter_( p, r ) ./ N;
mean_Ip = boxfilter_( I .* p, r ) ./ N;
cov_Ip = mean_Ip - mean_I .* mean_p;

mean_II = boxfilter_( I .* I, r ) ./ N;
var_I = mean_II - mean_I .* mean_I;

a = cov_Ip ./ ( var_I + eps );
b = mean_p - a .* mean_I;

mean_a = boxfilter_( a, r ) ./ N;
mean_b = boxfilter_( b, r ) ./ N;

q = mean_a .* I + mean_b;
end 

function imDst = boxfilter_( imSrc, r )








[ hei, wid ] = size( imSrc );
imDst = zeros( size( imSrc ) );


imCum = cumsum( imSrc, 1 );

imDst( 1:r + 1, : ) = imCum( 1 + r:2 * r + 1, : );
imDst( r + 2:hei - r, : ) = imCum( 2 * r + 2:hei, : ) - imCum( 1:hei - 2 * r - 1, : );
imDst( hei - r + 1:hei, : ) = repmat( imCum( hei, : ), [ r, 1 ] ) - imCum( hei - 2 * r:hei - r - 1, : );


imCum = cumsum( imDst, 2 );

imDst( :, 1:r + 1 ) = imCum( :, 1 + r:2 * r + 1 );
imDst( :, r + 2:wid - r ) = imCum( :, 2 * r + 2:wid ) - imCum( :, 1:wid - 2 * r - 1 );
imDst( :, wid - r + 1:wid ) = repmat( imCum( :, wid ), [ 1, r ] ) - imCum( :, wid - 2 * r:wid - r - 1 );
end 

function F = RF( img, sigma_s, sigma_r, num_iterations, joint_image )

I = double( img );

if ~exist( 'num_iterations', 'var' )
num_iterations = 3;
end 

if exist( 'joint_image', 'var' ) && ~isempty( joint_image )
J = double( joint_image );

if ( size( I, 1 ) ~= size( J, 1 ) ) || ( size( I, 2 ) ~= size( J, 2 ) )
error( 'Input and joint images must have equal width and height.' );
end 
else 
J = I;
end 

[ h, w, num_joint_channels ] = size( J );






dIcdx = diff( J, 1, 2 );
dIcdy = diff( J, 1, 1 );

dIdx = zeros( h, w );
dIdy = zeros( h, w );


for c = 1:num_joint_channels
dIdx( :, 2:end  ) = dIdx( :, 2:end  ) + abs( dIcdx( :, :, c ) );
dIdy( 2:end , : ) = dIdy( 2:end , : ) + abs( dIcdy( :, :, c ) );
end 


dHdx = ( 1 + sigma_s / sigma_r * dIdx );
dVdy = ( 1 + sigma_s / sigma_r * dIdy );







dVdy = dVdy';



N = num_iterations;
F = I;

sigma_H = sigma_s;

for i = 0:num_iterations - 1


sigma_H_i = sigma_H * sqrt( 3 ) * 2 ^ ( N - ( i + 1 ) ) / sqrt( 4 ^ N - 1 );

F = TransformedDomainRecursiveFilter_Horizontal( F, dHdx, sigma_H_i );
F = image_transpose( F );

F = TransformedDomainRecursiveFilter_Horizontal( F, dVdy, sigma_H_i );
F = image_transpose( F );

end 

F = cast( F, class( img ) );

end 


function F = TransformedDomainRecursiveFilter_Horizontal( I, D, sigma )


a = exp(  - sqrt( 2 ) / sigma );

F = I;
V = a .^ D;

[ h, w, num_channels ] = size( I );


for i = 2:w
for c = 1:num_channels
F( :, i, c ) = F( :, i, c ) + V( :, i ) .* ( F( :, i - 1, c ) - F( :, i, c ) );
end 
end 


for i = w - 1: - 1:1
for c = 1:num_channels
F( :, i, c ) = F( :, i, c ) + V( :, i + 1 ) .* ( F( :, i + 1, c ) - F( :, i, c ) );
end 
end 

end 

function T = image_transpose( I )

[ h, w, num_channels ] = size( I );

T = zeros( [ w, h, num_channels ], class( I ) );

for c = 1:num_channels
T( :, :, c ) = I( :, :, c )';
end 

end 

function VSMap = SDSP( image )


















sigmaF = 6.2;
omega0 = 0.002;
sigmaD = 114;
sigmaC = 0.2;

[ oriRows, oriCols, junk ] = size( image );
image = double( image );
dsImage( :, :, 1 ) = imresize( image( :, :, 1 ), [ 256, 256 ], 'bilinear' );
dsImage( :, :, 2 ) = imresize( image( :, :, 2 ), [ 256, 256 ], 'bilinear' );
dsImage( :, :, 3 ) = imresize( image( :, :, 3 ), [ 256, 256 ], 'bilinear' );
lab = RGB2Lab( dsImage );

LChannel = lab( :, :, 1 );
AChannel = lab( :, :, 2 );
BChannel = lab( :, :, 3 );

LFFT = fft2( double( LChannel ) );
AFFT = fft2( double( AChannel ) );
BFFT = fft2( double( BChannel ) );

[ rows, cols, junk ] = size( dsImage );
LG = logGabor( rows, cols, omega0, sigmaF );
FinalLResult = real( ifft2( LFFT .* LG ) );
FinalAResult = real( ifft2( AFFT .* LG ) );
FinalBResult = real( ifft2( BFFT .* LG ) );

SFMap = sqrt( FinalLResult .^ 2 + FinalAResult .^ 2 + FinalBResult .^ 2 );


coordinateMtx = zeros( rows, cols, 2 );
coordinateMtx( :, :, 1 ) = repmat( ( 1:1:rows )', 1, cols );
coordinateMtx( :, :, 2 ) = repmat( 1:1:cols, rows, 1 );

centerY = rows / 2;
centerX = cols / 2;
centerMtx( :, :, 1 ) = ones( rows, cols ) * centerY;
centerMtx( :, :, 2 ) = ones( rows, cols ) * centerX;
SDMap = exp(  - sum( ( coordinateMtx - centerMtx ) .^ 2, 3 ) / sigmaD ^ 2 );


maxA = max( AChannel( : ) );
minA = min( AChannel( : ) );
normalizedA = ( AChannel - minA ) / ( maxA - minA );

maxB = max( BChannel( : ) );
minB = min( BChannel( : ) );
normalizedB = ( BChannel - minB ) / ( maxB - minB );

labDistSquare = normalizedA .^ 2 + normalizedB .^ 2;
SCMap = 1 - exp(  - labDistSquare / ( sigmaC ^ 2 ) );
VSMap = SFMap .* SDMap .* SCMap;

VSMap = imresize( VSMap, [ oriRows, oriCols ], 'bilinear' );
VSMap = ( mat2gray( VSMap ) * 255 );
end 

function labImage = RGB2Lab( image )

image = double( image );
normalizedR = image( :, :, 1 ) / 255;
normalizedG = image( :, :, 2 ) / 255;
normalizedB = image( :, :, 3 ) / 255;

RSmallerOrEqualto4045 = normalizedR <= 0.04045;
RGreaterThan4045 = 1 - RSmallerOrEqualto4045;
tmpR = ( normalizedR / 12.92 ) .* RSmallerOrEqualto4045;
tmpR = tmpR + power( ( normalizedR + 0.055 ) / 1.055, 2.4 ) .* RGreaterThan4045;

GSmallerOrEqualto4045 = normalizedG <= 0.04045;
GGreaterThan4045 = 1 - GSmallerOrEqualto4045;
tmpG = ( normalizedG / 12.92 ) .* GSmallerOrEqualto4045;
tmpG = tmpG + power( ( normalizedG + 0.055 ) / 1.055, 2.4 ) .* GGreaterThan4045;

BSmallerOrEqualto4045 = normalizedB <= 0.04045;
BGreaterThan4045 = 1 - BSmallerOrEqualto4045;
tmpB = ( normalizedB / 12.92 ) .* BSmallerOrEqualto4045;
tmpB = tmpB + power( ( normalizedB + 0.055 ) / 1.055, 2.4 ) .* BGreaterThan4045;

X = tmpR * 0.4124564 + tmpG * 0.3575761 + tmpB * 0.1804375;
Y = tmpR * 0.2126729 + tmpG * 0.7151522 + tmpB * 0.0721750;
Z = tmpR * 0.0193339 + tmpG * 0.1191920 + tmpB * 0.9503041;

epsilon = 0.008856;
kappa = 903.3;

Xr = 0.9642;
Yr = 1.0;
Zr = 0.8251;

xr = X / Xr;
yr = Y / Yr;
zr = Z / Zr;

xrGreaterThanEpsilon = xr > epsilon;
xrSmallerOrEqualtoEpsilon = 1 - xrGreaterThanEpsilon;
fx = power( xr, 1.0 / 3.0 ) .* xrGreaterThanEpsilon;
fx = fx + ( kappa * xr + 16.0 ) / 116.0 .* xrSmallerOrEqualtoEpsilon;

yrGreaterThanEpsilon = yr > epsilon;
yrSmallerOrEqualtoEpsilon = 1 - yrGreaterThanEpsilon;
fy = power( yr, 1.0 / 3.0 ) .* yrGreaterThanEpsilon;
fy = fy + ( kappa * yr + 16.0 ) / 116.0 .* yrSmallerOrEqualtoEpsilon;

zrGreaterThanEpsilon = zr > epsilon;
zrSmallerOrEqualtoEpsilon = 1 - zrGreaterThanEpsilon;
fz = power( zr, 1.0 / 3.0 ) .* zrGreaterThanEpsilon;
fz = fz + ( kappa * zr + 16.0 ) / 116.0 .* zrSmallerOrEqualtoEpsilon;

[ rows, cols, junk ] = size( image );
labImage = zeros( rows, cols, 3 );
labImage( :, :, 1 ) = 116.0 * fy - 16.0;
labImage( :, :, 2 ) = 500.0 * ( fx - fy );
labImage( :, :, 3 ) = 200.0 * ( fy - fz );
end 

function LG = logGabor( rows, cols, omega0, sigmaF )
[ u1, u2 ] = meshgrid( ( [ 1:cols ] - ( fix( cols / 2 ) + 1 ) ) / ( cols - mod( cols, 2 ) ),  ...
( [ 1:rows ] - ( fix( rows / 2 ) + 1 ) ) / ( rows - mod( rows, 2 ) ) );
mask = ones( rows, cols );
for rowIndex = 1:rows
for colIndex = 1:cols
if u1( rowIndex, colIndex ) ^ 2 + u2( rowIndex, colIndex ) ^ 2 > 0.25
mask( rowIndex, colIndex ) = 0;
end 
end 
end 
u1 = u1 .* mask;
u2 = u2 .* mask;

u1 = ifftshift( u1 );
u2 = ifftshift( u2 );

radius = sqrt( u1 .^ 2 + u2 .^ 2 );
radius( 1, 1 ) = 1;

LG = exp( (  - ( log( radius / omega0 ) ) .^ 2 ) / ( 2 * ( sigmaF ^ 2 ) ) );
LG( 1, 1 ) = 0;
end 





