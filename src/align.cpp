#include "align.h"
#include <string>
#include <vector>
#include <cmath>
//grayscale!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
using std::string;
using std::cout;
using std::cerr;
using std::tie;
using std::endl;
using std::make_tuple;
Image gray_world(Image );
Matrix<double> kernel;
const int SHBOUNDS = 15;
const int eps = 7;
const int off = 15;
class BoxFilterOp
{	
  public:
    std::tuple<int, int, int> operator () (const Image &im) const{
		radius = kernel.n_rows/2;
		uint size = 2 * radius + 1;
        uint r, g, b;
		double norm = 0;
		double sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = im(i, j);
				sum_r += r*kernel(i,j);
                sum_g += g*kernel(i,j);
                sum_b += b*kernel(i,j);
            }
        }
//	auto norm = size*size;
    return make_tuple(int(sum_r), int(sum_g),int(sum_b));
    }
	static int radius;
    // Radius of neighbourhoud, which is passed to that operator
};
class MedianFilter
{	
  public:
    std::tuple<int, int, int> operator () (const Image &im) const{
		uint size = 2 * radius + 1;
        uint r, g, b;
		double norm = 0;
		std::vector<int> red,green,blue;	
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = im(i, j);
				red.push_back(r);
				green.push_back(g);
				blue.push_back(b); }
        }
		std::sort(red.begin(),red.end());
		std::sort(green.begin(),green.end());
		std::sort(blue.begin(),blue.end());

//	auto norm = size*size;
    return make_tuple(red[red.size()/2+1],green[green.size()/2+1],blue[blue.size()/2+1]);
    }
	static int radius;
    // Radius of neighbourhoud, which is passed to that operator
};
int BoxFilterOp::radius;
int MedianFilter::radius;
Image overlay2(const Image &Red, const Image &Green, const Image &Blue)
{
    Image ans(std::min(Red.n_rows,std::min(Green.n_rows,Blue.n_rows)),std::min(Red.n_cols,std::min(Green.n_cols,Blue.n_cols)));
    int i, j;
    for (i = 0; i < ans.n_rows; i++)
    {
		for (j = 0; j < ans.n_cols; j++)
		{
			if ((get<0>(Red(i,j))!=0) && (get<0>(Green(i,j))!=0) && (get<0>(Blue(i,j))!=0))
			    ans(i, j) = make_tuple(get<0>(Red(i, j)), get<0>(Green(i, j)), get<0>(Blue(i, j)));
		}
    }
    return ans;
}
Image shift(const Image & one,int n,int k){
		Image ans(one.n_rows,one.n_cols);//возможно обнулим
		n = -n;
		if (n<0)
			if (k>0){
				for (int j=0;(j < ans.n_rows) && (j<ans.n_rows+n);j++)
					for (int i=k;(i<ans.n_cols) && (i<ans.n_cols+k);i++){
						ans(j,i)=one(j-n,i-k);
					}
			}
			else{
				for (int j=0;(j < ans.n_rows) && (j<ans.n_rows+n);j++)
					for (int i=0;(i<ans.n_cols) && (i<ans.n_cols+k);i++){
						ans(j,i)=one(j-n,i-k);
					}
			}
		else{
			if (k>0){
				for (int j=n;(j < one.n_rows) && (j<one.n_rows+n);j++)
					for (int i=k;(i<ans.n_cols) && (i<ans.n_cols+k);i++){
						ans(j,i)=one(j-n,i-k);
					}
			}
			else{
				for (int j=n;(j < one.n_rows) && (j<one.n_rows+n);j++)
					for (int i=0;(i<ans.n_cols) && (i<ans.n_cols+k);i++){
						ans(j,i)=one(j-n,i-k);
					}
			}
		}
		return ans;
}


double MSE(const Image & one, const Image & two,int n,int k){
	double ans2=0;
		if (n<0)
			if (k>0){
				for (int j=0;(j<one.n_rows) && (j+k<two.n_rows);j++)
					for (int i=0;(i<one.n_cols) && (i-n<two.n_cols);i++){
						ans2+=((int)get<0>(one(j,i))-(int)get<0>(two(j+k,i-n)))*((int)get<0>(one(j,i))-(int)get<0>(two(j+k,i-n)));
					}
				 	return ans2/((two.n_cols+n)*(two.n_rows-k));
			}
			else{
				for (int j=0;(j-k<one.n_rows) && (j<two.n_rows);j++)
					for (int i=0;(i<one.n_cols) && (i-n<two.n_cols);i++){
						ans2+=((int)get<0>(one(j-k,i))-(int)get<0>(two(j,i-n)))*((int)get<0>(one(j-k,i))-(int)get<0>(two(j,i-n)));
					}
				return ans2/((one.n_rows+k)*(two.n_cols+n));
			}
		else{
			if (k>0){
				for (int j=0;(j<one.n_rows) && (j+k<two.n_rows);j++)
					for (int i = 0;(i+n<one.n_cols) && (i<two.n_cols);i++){
						ans2+=((int)get<0>(one(j,i+n))-(int)get<0>(two(j+k,i)))*((int)get<0>(one(j,i+n))-(int)get<0>(two(j+k,i)));
					}
			 	return ans2/((one.n_cols-n)*(two.n_rows-k));
			}
			else{ 
				for (int j=0;(j-k<one.n_rows) && (j<two.n_rows);j++)
					for (int i = 0;(i+n<one.n_cols) && (i<two.n_cols);i++){
						ans2+=((int)get<0>(one(j-k,i+n))-(int)get<0>(two(j,i)))*((int)get<0>(one(j-k,i+n))-(int)get<0>(two(j,i)));
					}
		 	return ans2/((one.n_cols-n)*(one.n_rows+k));					
			}
		}
}
std::tuple <int,int> minindex(const Matrix<double> & a){
	float minv=a(0,0);
	uint resi=0;
	uint resj=0;
	for (uint i=0;i<a.n_rows;++i){
		for (uint j=0;j<a.n_cols;++j){
			if (a(i,j)<minv){
				minv=a(i,j);
				resi=i;
				resj=j;
			}
		}
	}
	return make_tuple(resi-SHBOUNDS,resj-SHBOUNDS);
}
std::tuple<int,int> bestmetric(const Image & one, const Image & two){
	Matrix<double> metrics(2*SHBOUNDS+1,2*SHBOUNDS+1);

	for (int i=-SHBOUNDS;i<SHBOUNDS+1;i++){
		for (int j = -SHBOUNDS; j<SHBOUNDS+1;j++){
			metrics(j+SHBOUNDS,i+SHBOUNDS) = MSE(one,two,i,j);
		}
	}
	return minindex(metrics);    
}
int abs(int i){
	if (i<0)
		return -i;
	return i;
}
//n rows
Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
			bool isInterp, bool isSubpixel, double subScale)
{
	int i, j;
	if (isPostprocessing){
		if (postprocessingType == "--gray-world")
			srcImage = gray_world(srcImage);
	} 
//	srcImage = canny(srcImage, 70 , 100);
	int delim =srcImage.n_rows/3; 
	Image Blue(delim-off*2,srcImage.n_cols-off*2),Green(delim-off*2,srcImage.n_cols-off*2),Red(delim-off*2,srcImage.n_cols-off*2);
	for (i = 0; i < delim-off*2; i++) {
		for (j = off; j < srcImage.n_cols-off; j++) {
			Blue(i,j-off)=srcImage(i+off,j);
			Green(i,j-off)=srcImage(i+delim+off,j);
			Red(i,j-off)=srcImage(i+2*delim+off,j);
		}
	}
	if (isSubpixel){
			Blue = resize(Blue,subScale);
			Red = resize(Red,subScale);
			Green = resize(Green,subScale);		
			std::tie(i,j)=bestmetric(Green,Blue);
			Blue = shift(Blue,i,j);	
			std::tie(i,j)=bestmetric(Green,Red);	
			Red = shift(Red,i,j);
		return resize(overlay2(Red,Green,Blue),1/subScale);
//			Blue = resize(Blue,1/subScale);
//			Red = resize(Red,1/subScale);
//			Green = resize(Green,1/subScale);
//		return 	overlay2(Red,Green,Blue);
	}
	else{
	std::tie(i,j)=bestmetric(Green,Blue);

	Blue = shift(Blue,i,j);	
	std::tie(i,j)=bestmetric(Green,Red);
	Red = shift(Red,i,j);			
	}
	if (isPostprocessing){
		if (postprocessingType == "--unsharp")
			return unsharp(overlay2(Red,Green,Blue));
	}
	return overlay2(Red,Green,Blue);
}
Image sobel_x(Image src_image) {
	BoxFilterOp::radius = 1;
	kernel = {{-1, 0, 1},
							 {-2, 0, 2},
							 {-1, 0, 1}};

	return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
	BoxFilterOp::radius = 1;
	kernel = {{ 1,  2,  1},
							 { 0,  0,  0},
							 {-1, -2, -1}};
	return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
	kernel = {{-1.0/6, -2.0/3, -1.0/6},
			{-2.0/3,13.0/3,-2.0/3},
			{-1.0/6,-2.0/3,-1.0/6}};
	// kernel = {{-1.0/10,-2.0/10,-1.0/10},
	// 			{-2.0/10,22.0/10,-2.0/20},
	// 			{-1.0/10,-2.0/10,-1.0/10}};
	BoxFilterOp::radius = 1;
	src_image = custom(src_image,kernel);
	for (uint i = 0; i < src_image.n_rows; i++){
		for (uint j = 0; j < src_image.n_cols; j++){
			if (get<0>(src_image(i,j)) < 0)
				get<0>(src_image(i,j)) = 0;
			if (get<1>(src_image(i,j)) < 0)
				get<1>(src_image(i,j)) = 0;
			if (get<2>(src_image(i,j)) < 0)
				get<2>(src_image(i,j)) = 0;
			if (get<0>(src_image(i,j)) > 255)
				get<0>(src_image(i,j)) = 255;
			if (get<1>(src_image(i,j)) > 255)
				get<1>(src_image(i,j)) = 255;
			if (get<2>(src_image(i,j)) > 255)
				get<2>(src_image(i,j)) = 255;				
		}
	}
	return src_image;
}
double mean(const Image & im){
	double ans=0;
	for (uint i = 0;i < im.n_rows; i++){
		for (uint j=0;j < im.n_cols;j++){
			ans+=get<0>(im(i,j));
		}
	}
	return ans/(im.n_cols*im.n_rows);
}
Image gray_world(Image src_image) {
	double sr=0,sg=0,sb=0;
	uint i,j;

	double delim = src_image.n_rows/3;
	for (i = off ; i < src_image.n_rows/3 - off; i++) {
		for (j = off; j < src_image.n_cols - off; j++) {
			sb+=get<0>(src_image(i,j));	
		}
	}
		sb=sb/((delim-2*off)*(src_image.n_cols-2*off));	
	for (i = src_image.n_rows/3 + off; i < src_image.n_rows*2/3 - off; i++) {
		for (j = off; j < src_image.n_cols-off; j++) {
			sg+=get<0>(src_image(i,j));
		}
	}
		sg=sg/((delim-2*off)*(src_image.n_cols-2*off));
	for (i = src_image.n_rows*2/3+off; i < src_image.n_rows-off; i++) {
		for (j = off; j < src_image.n_cols-off; j++) {
			sr+=get<0>(src_image(i,j));
		}
	}
		sr=sr/((delim-2*off)*(src_image.n_cols-2*off));
	double s = (sr+sg+sb)/3;
	sr=s/sr;
	sb=s/sb;
	sg=s/sg;

	for (i = off; i < src_image.n_rows/3-off; i++) {
		for (j = off; j < src_image.n_cols-off; j++) {
			uint koef = double(get<0>(src_image(i,j)))*sb;
			if (koef > 255)
				koef = 255;
			src_image(i,j)=make_tuple(koef,koef,koef);
		}
	}	
	for (i = src_image.n_rows/3+off; i < src_image.n_rows*2/3-off; i++) {
		for (j = off; j < src_image.n_cols-off; j++) {
			uint koef = double(get<0>(src_image(i,j)))*sg;
			if (koef > 255)
				koef = 255;
			src_image(i,j)=make_tuple(koef,koef,koef);		}
	}
	for (i = src_image.n_rows*2/3 + off; i < src_image.n_rows-off; i++) {
		for (j = off; j < src_image.n_cols-off; j++) {
			uint koef = double(get<0>(src_image(i,j)))*sr;
			if (koef > 255)
				koef = 255;
			src_image(i,j)=make_tuple(koef,koef,koef);
		}
	}
	return src_image;
}

Image resize(Image src_image, double scale) {
	uint i,j,k,l;
	if (scale > 1){
		Image Big(src_image.n_rows*scale,src_image.n_cols*scale);
		for (i = 0; i < src_image.n_rows; i++){
			for (j = 0; j < src_image.n_cols; j++){
				for (k = 0; k < scale; k++){
					for (l=0 ; l < scale ;l++)
						Big(i*scale+k,j*scale+l) = src_image(i,j);  	
				}
			}
		}	
		return Big;
	}

	Image Little(double(src_image.n_rows)*scale,double(src_image.n_cols)*scale);
	for (i = 0; i < int(src_image.n_rows); i += (1.0/scale))
		for (j = 0; j < int(src_image.n_cols); j += (1.0/scale))
			Little(i*scale,j*scale) = src_image(i,j);
	return Little;
}

Image custom(Image src_image, Matrix<double> kernel1) {
	// Function custom is useful for making concrete linear filtrations
	// like gaussian or sobel. So, we assume that you implement custom
	// and then implement other filtrations using this function.
	// sobel_x and sobel_y are given as an example.
	return src_image.unary_map(BoxFilterOp());
//	return src_image;

}

Image autocontrast(Image src_image, double fraction) {
	return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
	BoxFilterOp::radius = radius;
	kernel = Matrix<double>(2*radius+1,2*radius+1);
	int x,y ;
	for (x = -radius; x < radius+1; x++){
		for (y = -radius; y < radius+1;y++){
			kernel(x+radius,y+radius)=exp(-(x*x+y*y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
		}
	}
	return custom(src_image,kernel);
}
Image gaussian_separable(Image src_image, double sigma, int radius) {
	int k,x,y ;
	uint i,j;
	std::vector<double> ver(2*radius+1);
	Matrix<double> ab(2*radius+1,2*radius+1);
	for (x = -radius; x < radius+1; x++){
		for (y = -radius; y < radius+1;y++){
			ab(x+radius,y+radius) = exp(-(x*x+y*y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
		}
	}

	for (x = 0;x < 2*radius+1; x++){
		ver[x] = 0.0;
		for(y=0; y < 2*radius + 1;y++){
			ver[x] += ab(x,y); 	
		}

	}
	Image out(src_image.n_rows,src_image.n_cols);

	double sumr = 0,sumb = 0,sumg = 0;
	for (i=0;i<src_image.n_rows;i++){
		for (j=radius;j<src_image.n_cols-radius;j++){
			for ( k = -radius ;k < radius ; k++){
				sumr += get<0>(src_image(i,j+k))*ver[k+radius]; 
				sumg += get<1>(src_image(i,j+k))*ver[k+radius];
				sumb += get<2>(src_image(i,j+k))*ver[k+radius];
			}
		out(i,j) = make_tuple(sumr,sumg,sumb);
			sumr = 0;
			sumb = 0;
			sumg = 0;
		}
	}
	sumr = 0;sumb=0;sumg=0;
	src_image = out.deep_copy();
	for (j=0;j<src_image.n_cols;j++){
		for (i=radius;i<src_image.n_rows-radius;i++){
			for ( k = -radius ;k < radius ; k++){
				sumr += get<0>(src_image(i+k,j))*ver[k+radius]; 
				sumg += get<1>(src_image(i+k,j))*ver[k+radius];
				sumb += get<2>(src_image(i+k,j))*ver[k+radius];
			}
		out(i,j) = make_tuple(sumr,sumg,sumb);
			sumr= 0;
			sumb = 0;
			sumg = 0;
		}
	}
	return out;
}
Image median(Image src_image, int radius) {
	MedianFilter::radius = radius;	
	return src_image.unary_map(MedianFilter());
}

Image median_linear(Image src_image, int radius) {
	return src_image;
}

Image median_const(Image src_image, int radius) {
	return src_image;
}
void Fill(Matrix<double> & labels, int x, int y,Matrix<int> & used)
{
	if ((x < 0) || (x == labels.n_rows))
		return;
	if ((y < 0) || (y == labels.n_cols))
		return;
	if (!labels(x,y))
		return;
	used(x,y) = 1;
//	cout << "xy" << x << " " << y << endl;
	if (labels(x,y) == 127)
		labels(x,y) = 255;
	if (y > 0 )
	if	(!used(x,y-1))
		Fill(labels, x, y - 1, used);
	if (y < labels.n_cols-1)
		if (!used(x,y+1))
		Fill(labels, x, y + 1,  used);
	if (x > 0 )
		if (!used(x-1,y))
		Fill(labels, x - 1, y, used);
	if (x < labels.n_rows-1)
		if (!used(x+1,y))
		Fill(labels, x + 1, y, used);
	if ((x < labels.n_rows-1) && (y < labels.n_cols-1))
		if (!used(x+1,y+1))
		Fill(labels, x + 1 , y + 1 , used);
	if ((x > 0 ) && (y > 0))
		if (!used(x-1,y-1))
		Fill(labels, x - 1 , y - 1 , used);
	if ((x > 0 ) && (y < labels.n_cols-1))
		if (!used(x-1,y+1))
		Fill(labels, x - 1 , y + 1 , used);
	if ((x < labels.n_rows-1) && (y >0))
		if (!used(x+1,y-1))
		Fill(labels, x + 1 , y - 1 , used);
}
void Labeling(Matrix<double> & labels,Matrix<int> & used){
	for (int x = 0; x < labels.n_rows; x++)
		for (int y = 0; y < labels.n_cols; y++)
		{
			if (labels(x, y) == 255)
				Fill(labels, x, y,used);
		}
}
std::tuple<uint,uint> Findindex(Matrix<double> map,uint starti,uint endi,uint startj,uint endj){
	double sum=0;
	double max=-1;
	int maxi=0;
	uint i,j;
	for (i = starti; i < endi; i++ ){
		for( j = startj ; j < endj; j++){
			sum += map(i,j);
		}
		if (sum > max){
			max = sum;
			maxi = i;
		}

		sum=0;
	}

	if (maxi - eps <0)
	for (i = 0; (i <= maxi + eps) && (i < map.n_rows) ; i++){
			for( j = startj ; j < endj; j++){
				map(i,j)=0;
			}

	}
	else
	for (i=maxi-eps;(i <= maxi + eps) && (i < map.n_rows) ;i++){
			for( j = startj ; j < endj; j++){
				map(i,j)=0;
			}

	}
	int newmaxi = 0;
	sum = 0;
	max = -1;
	for (i = starti; i < endi; i++ ){
			for( j = startj; j < endj; j++){
				sum += map(i,j);
			}
			if (sum > max){
				max = sum;
				newmaxi = i;
			}
			sum=0;
	}
	return make_tuple(maxi,newmaxi);
}
std::tuple<uint,uint> Findindexh(Matrix<double> map,uint starti,uint endi,uint startj,uint endj){
	double sum=0;
	double max=0;
	int maxi=0;
	uint i,j;
	for( j = startj ; j < endj; j++){
		for (i = starti; i < endi; i++ ){
			sum += map(i,j);
		}
		if (sum > max){
			max = sum;
			maxi = j;
		}
		sum = 0;
	}
	if (maxi - eps <0)
	for( j = 0 ; j <= maxi+eps && j < map.n_cols; j++){
		if (j!=maxi){
			for (i=starti;i<endi;i++){
				map(i,j)=0;
			}
		}
	}
	else
	for( j = maxi-eps ; j <= maxi+eps && j < map.n_cols; j++){
		if (j!=maxi){
			for (i=starti;i<endi;i++){
				map(i,j)=0;
			}
		}
	}
	int newmaxi = 0;
	max = 0;
	sum = 0;
	for( j = startj; j < endj; j++){
		if (j!=maxi){
		for (i = starti; i < endi; i++ ){
				sum += map(i,j);
			}
			if (sum > max){
				max = sum;
				newmaxi = j;
			}
			sum = 0;	
		}
	}
	return make_tuple(maxi,newmaxi);
}
Image canny(Image src_image, int threshold1, int threshold2) {
	Image zz(src_image);
	zz = gaussian_separable(src_image,1.4,2);
	
	Image pox = sobel_x(zz);
	Image poy = sobel_y(zz);
	Matrix<double> map(src_image.n_rows,src_image.n_cols);
	uint i,j;

	for (i = 0; i < src_image.n_rows; i++){
		for (j = 0; j < src_image.n_cols; j++){
			map(i,j) = sqrt(get<0>(pox(i,j))*get<0>(pox(i,j))+get<0>(poy(i,j))*get<0>(poy(i,j)));
		}
	}	

	for (i=1; i < src_image.n_rows - 1; i++){
		for (j=1; j < src_image.n_cols - 1;j++){
			int alpha = round(atan2(get<0>(pox(i,j)),get<0>(poy(i,j)))*4/M_PI);
			if (alpha == 0 || alpha == 4 || alpha == -4){
				if (map(i,j) < map(i,j+1) || map(i,j) < map(i,j-1))
					map(i,j) = 0;
			}
			else if (alpha == 1 || alpha == -3){
				if (map(i,j) < map(i+1,j+1) || map(i,j) < map(i-1,j-1))
					map(i,j) = 0;			
			}
			else if (alpha == 2 || alpha == -2)
				if (map(i,j) < map(i+1,j) || map(i,j) < map(i-1,j))
					map(i,j) = 0;								
			else if (alpha == -1 || alpha == 3)
				if (map(i,j) < map(i+1,j-1) || map(i,j) < map(i-1,j+1))
					map(i,j) = 0;								
		}
	}
	const int high = 255,low = 127,clear = 0;
	for (i = 0; i < src_image.n_rows; i++){
		for (j = 0; j < src_image.n_cols; j++){
			if (map(i,j)<threshold1)
				map(i,j) = clear;
			else if (map(i,j)> threshold2)
				map(i,j) = high;
			else map(i,j) = low;		
		}
	}
//	map = vid(map);
	Matrix<int> used(map.n_rows,map.n_cols);
	for( uint i = 0; i < used.n_rows; i++ )
		for( uint j = 0; j < used.n_cols; j++ )
			used(i,j) = false;
	Labeling(map,used);

	for (i = 0; i < map.n_rows; i++){
	 	for (j = 0; j < map.n_cols; j++){
			 if (map(i,j)==127)
			 	map(i,j) = 0;
		 }
	}
	
	for (i = 0; i < src_image.n_rows; i++){
		for (j = 0; j < src_image.n_cols; j++){	
			src_image(i,j) = make_tuple(map(i,j),map(i,j),map(i,j)); 
		}
	}

	const double percent = 0.05;
	uint h = int(percent*double(map.n_rows));
	uint starti,endi,startj,endj;
	tie(i,j) = Findindex(map,3,h,0,map.n_cols);
	starti = std::max(i,j);
	tie(i,j) = Findindex(map,map.n_rows - h,map.n_rows-3,0,map.n_cols);
	endi = std::min(i,j);
	tie(i,j) = Findindexh(map,0,map.n_rows,3,h);
	startj = std::max(i,j);
	tie(i,j) = Findindexh(map,0,map.n_rows,map.n_cols-h,map.n_cols-3);
	endj =  std::min(i,j);
//	return src_image.submatrix(starti,startj,endi-starti,endj-startj);
	return src_image;
}


