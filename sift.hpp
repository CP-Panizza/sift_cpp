#ifndef _SIFT_HPP_
#define _SIFT_HPP_

#include <stdexcept>
#include <fstream>
#include <cmath>
#include <c++/cfloat>
#include "matrix.hpp"
#include "gaussfilter.hpp"

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct {
    int height;
    int width;
    double **Data;
} Mat;


typedef struct {
    int row;
    int col;
} Point;


/*
 * 构建汉森矩阵
 */
Mat *BuildHessian(Point *point, Matrix *img) {
    Mat *HessianMat;
    const int size = 2;
    int i;
    double f0, f1, f2, f3, f4, f5, f6, f7, f8;


    HessianMat = new Mat;
    HessianMat->width = size;
    HessianMat->height = size;

    HessianMat->Data = new double *[size];
    for (i = 0; i < size; i++) {
        HessianMat->Data[i] = new double[size];
    }

    f0 = img->data[point->row][point->col];//不用考虑指针越界，因为留下了一的范围
    f1 = img->data[point->row][point->col + 1];
    f2 = img->data[point->row + 1][point->col];
    f3 = img->data[point->row][point->col - 1];
    f4 = img->data[point->row - 1][point->col];
    f5 = img->data[point->row - 1][point->col + 1];
    f6 = img->data[point->row + 1][point->col + 1];
    f7 = img->data[point->row + 1][point->col - 1];
    f8 = img->data[point->row - 1][point->col - 1];

    HessianMat->Data[0][0] = (f1 + f3 - 2 * f0);
    HessianMat->Data[1][1] = (f2 + f4 - 2 * f0);
    HessianMat->Data[0][1] = (f8 + f6 - f5 - f7) / 4;
    HessianMat->Data[1][0] = (f8 + f6 - f5 - f7) / 4;

    return HessianMat;
}


void DeleteHessian(Mat *val) {
    for (int i = 0; i < val->height; ++i) {
        delete[](val->Data[i]);
    }
    delete[](val->Data);
}


typedef struct {
    double norm;
    double direction;
} Gradient;


typedef struct {
    Point point;
    int octavenum;
    int intervalnum;
    double sigma;//尺度数,先用图片顺序记着
    Gradient grad;//主方向0到350°
    double *descriptor;
    bool flag;
} KeyPoint;


typedef struct {
    KeyPoint *point;
    int number;
} KeyPointGroup;

typedef struct {
    Point *Img1;
    Point *Img2;
    int Matchnumber;
} MatchPoint;


class ImagePyramid {
public:
    int OctaveNum;    //金字塔层数
    int IntervalNum;  //每一层矩阵数
    Matrix **mats;    //金字塔数据
    int n;            //n + 3

    ImagePyramid() {}

    /**
     * 建立高斯金字塔
     * @param origin
     * @param t
     * @param s
     */
    void build(Matrix *origin, int t, int s) {
        this->n = s;
        this->OctaveNum = static_cast<int>(
                (log(double(min(origin->height, origin->width))) - log(double(t))) / log(2.0) + 1);
        this->IntervalNum = s + 3;
        float sigma = 0.5;
        float k = static_cast<float>(pow((float) 2, 1 / (float) this->IntervalNum));
        float filterSize;
        this->mats = new Matrix *[this->OctaveNum];
        for (int i = 0; i < this->OctaveNum; ++i) {
            this->mats[i] = new Matrix[this->IntervalNum];
        }

        for (int j = 0; j < this->OctaveNum; ++j) {
            if (j == 0) {
                this->mats[j][0] = *origin->Copy();
            } else {
                this->mats[j][0] = *((this->mats[j - 1][this->IntervalNum - 3]).DownSampling());
            }
            filterSize = static_cast<float>(pow((double) 2, j) * sigma);
            for (int i = 1; i < this->IntervalNum; ++i) {
                auto *filter = new GaussFilter(filterSize);
                this->mats[j][i] = *(filter->Convolution(&(this->mats[j][0])));
                delete (filter);
                filterSize = filterSize * k;
            }
        }
    }

    /**
     * 建立高斯差分金字塔
     */
    ImagePyramid *buildDiff() {
        ImagePyramid *diff = new ImagePyramid;
        diff->n = this->n;
        diff->OctaveNum = this->OctaveNum;
        diff->IntervalNum = this->IntervalNum - 1; //张数减1
        diff->mats = new Matrix *[diff->OctaveNum];
        for (int i = 0; i < diff->OctaveNum; ++i) {
            diff->mats[i] = new Matrix[diff->IntervalNum];
        }

        for (int j = 0; j < diff->OctaveNum; ++j) {
            for (int i = 0; i < diff->IntervalNum; ++i) {
                diff->mats[j][i] = *this->difference(&this->mats[j][i + 1], &this->mats[j][i]);
            }
        }

        return diff;
    }


    /**
     * 计算差分金字塔同一层中相邻的两个矩阵的差值矩阵
     * @param a
     * @param b
     * @return
     */
    Matrix *difference(Matrix *a, Matrix *b) {
        Matrix *diff = new Matrix;
        diff->height = a->height;
        diff->width = a->width;

        diff->data = new double *[diff->height];
        for (int i = 0; i < diff->height; ++i) {
            diff->data[i] = new double[diff->width];
        }

        for (int j = 0; j < diff->height; ++j) {
            for (int i = 0; i < diff->width; ++i) {
                diff->data[j][i] = a->data[j][i] - b->data[j][i];
            }
        }

        return diff;
    }
};


Gradient *GetGradient(Point *point, Matrix *img) {
    Gradient *grad;
    double dx;
    double dy;
    double angle;
    double f1, f3, f2, f4;

    grad = new Gradient;

    f3 = img->data[point->row + 1][point->col];
    f1 = img->data[point->row - 1][point->col];
    f4 = img->data[point->row][point->col + 1];
    f2 = img->data[point->row][point->col - 1];

    dx = (f3 - f1) / 2.0;
    dy = (f4 - f2) / 2.0;
    angle = atan2(dy, dx) * 180 / 3.1415926; //梯度方位角
    if (angle < 0) {
        angle = angle + 360;
    }

    grad->norm = sqrt(pow(dx, 2) + pow(dy, 2));//梯度模值
    grad->direction = angle;

    return grad;
}


int MaxIndex(double *array, int length) {
    int i;
    double max;
    int max_sign;

    max = array[0];
    max_sign = 0;

    for (i = 1; i < length; i++) {
        if (array[i] > max) {
            max = array[i];
            max_sign = i;
        }
    }

    return max_sign;
}

void adjustLocalExtrema(ImagePyramid *DoG, int o, int s, int x, int y, double contrastThreshold, double edgeThreshold,
                        double sigma, int n, int SIFT_FIXPT_SCALE) {
    double SIFT_MAX_INTERP_STEPS = 5;
    double SIFT_IMG_BORDER = 5;
    double img_scale = 1.0 / (255 * SIFT_FIXPT_SCALE);
    double deriv_scale = img_scale * 0.5;
    double second_deriv_scale = img_scale;
    double cross_deriv_scale = img_scale * 0.25;
    auto img = DoG->mats[o][s];
    auto prev = DoG->mats[o][s - 1];
    auto next = DoG->mats[o][s + 1];
    int i = 0;
    while (i < SIFT_MAX_INTERP_STEPS) {
        if (s < 1 || s > n || y < SIFT_IMG_BORDER || y >= img.width - SIFT_IMG_BORDER || x < SIFT_IMG_BORDER ||
            x >= img.height - SIFT_IMG_BORDER) {
            return;
        }
        double dD[] = {(img.data[x][y + 1] - img.data[x][y - 1]) * deriv_scale,
                       (img.data[x + 1][y] - img.data[x - 1][y]) * deriv_scale,
                       (next.data[x][y] - prev.data[x][y]) * deriv_scale};

        double v2 = (img.data[x][y]) * 2;
        double dxx = (img.data[x][y + 1] + img.data[x][y - 1] - v2) * second_deriv_scale;
        double dyy = (img.data[x + 1][y] + img.data[x - 1][y] - v2) * second_deriv_scale;
        double dss = (next.data[x][y] + prev.data[x][y] - v2) * second_deriv_scale;
        double dxy =
                (img.data[x + 1][y + 1] - img.data[x + 1][y - 1] - img.data[x - 1][y + 1] + img.data[x - 1][y - 1]) *
                cross_deriv_scale;
        double dxs = (next.data[x][y + 1] - next.data[x][y - 1] - prev.data[x][y + 1] + prev.data[x][y - 1]) *
                     cross_deriv_scale;
        double dys = (next.data[x + 1][y] - next.data[x - 1][y] - prev.data[x + 1][y] + prev.data[x - 1][y]) *
                     cross_deriv_scale;

        double H[3][3] = {{dxx, dxy, dxs},
                          {dxy, dyy, dys},
                          {dxs, dys, dss}};



    }
}


/**
 * 寻找关键点
 * @param ip
 * @param diff_ip
 * @return
 */
KeyPointGroup *SearchKeyPointPosition(ImagePyramid *ip, ImagePyramid *diff_ip) {
    int big_then_count;//大于计数
    int small_then_count; //小于计数
    int p, q;
    int octavenum = diff_ip->OctaveNum;
    int intervalnum = diff_ip->IntervalNum;
    double threshold2 = 0.04 / diff_ip->n; //极值点低对比度极值 0.04/n
    const int SIFT_FIXPT_SCALE = 1;
    double prelim_contr_thr = 0.5 * 0.04 / (diff_ip->n * 255 * SIFT_FIXPT_SCALE);
    int length = 1;//极值点个数, 第一个点没用
    double k = pow((float) 2, 1 / (float) intervalnum);
    KeyPointGroup *Keypointgroup = (KeyPointGroup *) malloc(sizeof(KeyPointGroup));
    KeyPoint *TempPoint = (KeyPoint *) malloc(sizeof(KeyPoint));
    KeyPoint *TempPoint2 = (KeyPoint *) malloc(sizeof(KeyPoint));
    TempPoint[0].point.row = 0;//第一个点没用
    TempPoint[0].point.col = 0;
    TempPoint2[0].point.row = 0;//第一个点没用
    TempPoint2[0].point.col = 0;


    for (int i = 0; i < octavenum; ++i) {
        for (int j = 1; j < intervalnum - 1; ++j) {
            for (int h = 1; h < diff_ip->mats[i][0].height - 1; ++h) {
                for (int w = 1; w < diff_ip->mats[i][0].width - 1; ++w) {
                    if (fabs(diff_ip->mats[i][j].data[h][w]) < prelim_contr_thr) { //对像素点进行阈值化去除
                        continue;
                    }
                    /**
                     * 三个尺度内的26个点进行极值点比较
                     */
                    big_then_count = 0;
                    small_then_count = 0;
                    for (int m = -1; m < 2; ++m) {
                        for (p = -1; p < 2; p++) {
                            for (q = -1; q < 2; q++) {
                                if (diff_ip->mats[i][j].data[h][w] >= diff_ip->mats[i][j + m].data[h + p][w + q]) {
                                    big_then_count++;
                                } else if (diff_ip->mats[i][j].data[h][w] <=
                                           diff_ip->mats[i][j + m].data[h + p][w + q]) {
                                    small_then_count++;
                                }
                            }
                        }
                    }
                    //此点是最大值或者最小值
                    if (big_then_count == 27 || small_then_count == 27) {
                        length++;
                        TempPoint = (KeyPoint *) realloc(TempPoint, sizeof(KeyPoint) * length);
                        TempPoint[length - 1].point.row = h;
                        TempPoint[length - 1].point.col = w;
                        TempPoint[length - 1].octavenum = i;
                        TempPoint[length - 1].intervalnum = j;
                        TempPoint[length - 1].sigma = pow((float) 2, i + j / (float) intervalnum);//先用图片顺序记
                    }
                    big_then_count = 0;
                    small_then_count = 0;
                }
            }
        }
    }


    double threshold = pow((double) 10.0, 2.0) / 10.0;
    Mat *Hessian;
    double Tr2, Det;
    int length2 = 1;
    //消除边缘点
    for (int i = 1; i < length; i++) {
        Hessian = BuildHessian(&TempPoint[i].point, &(ip->mats[TempPoint[i].octavenum][TempPoint[i].intervalnum]));
        Tr2 = pow(Hessian->Data[0][0] + Hessian->Data[1][1], 2.0);
        Det = Hessian->Data[0][0] * Hessian->Data[1][1] - Hessian->Data[1][0] * Hessian->Data[1][0];
        if ((Tr2 / Det) < threshold) {
            length2++;
            TempPoint2 = (KeyPoint *) realloc(TempPoint2, sizeof(KeyPoint) * length2);
            TempPoint2[length2 - 1].point.row = TempPoint[i].point.row;
            TempPoint2[length2 - 1].point.col = TempPoint[i].point.col;
            TempPoint2[length2 - 1].octavenum = TempPoint[i].octavenum;;
            TempPoint2[length2 - 1].intervalnum = TempPoint[i].intervalnum;
            TempPoint2[length2 - 1].sigma = TempPoint[i].sigma;//先用图片顺序记着
        }
        DeleteHessian(Hessian);
    }

    double sigma = 0.5;
    double direction[36];
    double FilterSize;
    GaussFilter *gaussfilter;
    int height, width;
    int WindowSize;
    Point *temp_point = new Point;
    Gradient *gradient;
    int DirectionIndex;
    int PrincipalDirectionIndex;
    //寻找主方向
    for (int i = 1; i < length2; i++) {

        for (int j = 0; j < 36; j++) {
            direction[j] = 0;
        }

        FilterSize = pow((double) 2, (int) TempPoint2[length2 - 1].octavenum) * sigma *
                     pow((double) k, (int) TempPoint2[length2 - 1].intervalnum);
        gaussfilter = new GaussFilter(FilterSize);
        WindowSize = gaussfilter->Size;
        height = ip->mats[TempPoint[i].octavenum][TempPoint[i].intervalnum].height;
        width = ip->mats[TempPoint[i].octavenum][TempPoint[i].intervalnum].width;

        for (int m = -(WindowSize - 1) / 2; m < (WindowSize - 1) / 2; m++) {
            for (int n = -(WindowSize - 1) / 2; n < (WindowSize - 1) / 2; n++) {
                if (TempPoint2[i].point.row + m < 1
                    || TempPoint2[i].point.row + m > height - 2
                    || TempPoint2[i].point.col + n < 1
                    || TempPoint2[i].point.col + n > width - 2) {
                    continue;//如果超出图像范围以及求差分的范围，则跳出，不计算该点的梯度
                }
                temp_point->row = TempPoint2[i].point.row + m;
                temp_point->col = TempPoint2[i].point.col + n;
                gradient = GetGradient(temp_point, &(ip->mats[TempPoint[i].octavenum][TempPoint[i].intervalnum]));
                DirectionIndex = static_cast<int>(gradient->direction / 10);//这样描述不准确
                direction[DirectionIndex] = direction[DirectionIndex] + gradient->norm * gaussfilter->FilterData[m +
                                                                                                                 (WindowSize -
                                                                                                                  1) /
                                                                                                                 2][n +
                                                                                                                    (WindowSize -
                                                                                                                     1) /
                                                                                                                    2];
                delete (gradient);
            }
        }
        delete (gaussfilter);
        PrincipalDirectionIndex = MaxIndex(direction, 36);

        TempPoint2[i].grad.norm = direction[PrincipalDirectionIndex];
        TempPoint2[i].grad.direction = PrincipalDirectionIndex * 10;
    }

    Keypointgroup->point = TempPoint2;
    Keypointgroup->number = length2;

    return Keypointgroup;

}

typedef struct {
    int *angle;
    double *norm;
} GradHist;

typedef struct {
    Point LefttopPoint;
    int width;
    int height;
} Rect;

typedef struct {
    Matrix *NormImage;
    Matrix *DirectionImage;
} GradientImage;


Matrix *CropImage(Matrix *Img, Rect *Rect) {
    Matrix *CropImg;
    int i;
    int j;
    int k;
    int height;
    int width;
    int channels;
    int step;

    CropImg = new Matrix;

    height = Rect->height;
    width = Rect->width;
    CropImg->height = height;
    CropImg->width = width;

    CropImg->data = new double *[height];
    for (i = 0; i < height; i++) {
        CropImg->data[i] = new double[width];
    }


    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            CropImg->data[i][j] = Img->data[Rect->LefttopPoint.row + i][Rect->LefttopPoint.col + j];
        }
    }

    return CropImg;
}


double BilinearInterpolation(double pixel[4], Point p) {
    double value;

    value = pixel[0] * (1 - p.row) * (1 - p.col) + pixel[1] * p.col * (1 - p.row) + pixel[2] * p.row * (1 - p.col) +
            pixel[3] * p.row * p.col;

    return value;
}


Matrix *RotateImage(Matrix *Img, int MainDerection) {
    Matrix *RotateImg;
    Point p;
    int i;
    int j;
    int k;
    int height;
    int width;
    double x, y;
    double xr, yr;
    double pixel[4];
    double XAxisLength;
    double YAxisLength;

    RotateImg = new Matrix;

    height = Img->height;
    width = Img->width;
    RotateImg->height = height;
    RotateImg->width = width;
    XAxisLength = ((double) width - 1) / 2;
    YAxisLength = ((double) height - 1) / 2;

    RotateImg->data = new double *[height];
    for (i = 0; i < height; i++) {
        RotateImg->data[i] = new double[width];
    }

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            xr = j - XAxisLength;
            yr = -i + YAxisLength;
            x = cos((double) MainDerection * 3.1415926 / 180.0) * xr -
                sin((double) MainDerection * 3.1415926 / 180.0) * yr;
            y = sin((double) MainDerection * 3.1415926 / 180.0) * xr +
                cos((double) MainDerection * 3.1415926 / 180.0) * yr;
            //超出范围则赋值0
            if ((int) (x + XAxisLength) < 0 || (int) (-y + YAxisLength) < 0
                || (int) (x + XAxisLength) >= width - 1 || (int) (-y + YAxisLength) >= height - 1) {
                RotateImg->data[i][j] = 0;
                continue;
            }
            //双线性插值得到旋转图像
            p.col = x - (int) x;
            p.row = y - (int) y;
            pixel[0] = Img->data[(int) (-y + YAxisLength)][(int) (x + XAxisLength)];
            pixel[1] = Img->data[(int) (-y + YAxisLength)][(int) (x + XAxisLength + 1)];
            pixel[2] = Img->data[(int) (-y + YAxisLength + 1)][(int) (x + XAxisLength)];
            pixel[3] = Img->data[(int) (-y + YAxisLength + 1)][(int) (x + XAxisLength + 1)];
            RotateImg->data[i][j] = BilinearInterpolation(pixel, p);
        }
    }

    return RotateImg;

}


GradientImage *GetGradientImage(Matrix *Img) {
    GradientImage *GradientImg;
    Point p;
    Gradient *grad;
    int i;
    int j;
    int k;
    int height;
    int width;

    GradientImg = new GradientImage;
    GradientImg->NormImage = new Matrix;
    GradientImg->DirectionImage = new Matrix;

    height = Img->height;
    width = Img->width;
    GradientImg->NormImage->height = height;
    GradientImg->NormImage->width = width;
    GradientImg->DirectionImage->height = height;
    GradientImg->DirectionImage->width = width;


    GradientImg->NormImage->data = new double *[height];
    for (i = 0; i < height; i++) {
        GradientImg->NormImage->data[i] = new double[width];
    }

    GradientImg->DirectionImage->data = new double *[height];
    for (i = 0; i < height; i++) {
        GradientImg->DirectionImage->data[i] = new double[width];
    }

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            //边缘无法计算梯度，梯度置为0
            if (i - 1 < 0 || j - 1 < 0 || i + 2 > height || j + 2 > width) {
                GradientImg->DirectionImage->data[i][j] = 0;
                GradientImg->NormImage->data[i][j] = 0;
                continue;
            }
            p.row = i;
            p.col = j;
            grad = GetGradient(&p, Img);
            GradientImg->DirectionImage->data[i][j] = grad->direction;
            GradientImg->NormImage->data[i][j] = grad->norm;
            delete (grad);
        }
    }

    return GradientImg;
}


GradHist *GetGradHist(GradientImage *GradientImg, int AngularResolution, Rect *rect) {
    GradHist *Gradhist;
    int AngleNumber;
    int temp;
    int i, j;
    int index;

    AngleNumber = 360 / AngularResolution;

    Gradhist = (GradHist *) malloc(sizeof(GradHist));
    Gradhist->angle = (int *) malloc(sizeof(int) * AngleNumber);
    Gradhist->norm = (double *) malloc(sizeof(double) * AngleNumber);

    for (i = 0; i < AngleNumber; i++) {
        Gradhist->angle[i] = i * AngularResolution;
        Gradhist->norm[i] = 0;
    }

    for (i = rect->LefttopPoint.row; i < rect->LefttopPoint.row + rect->height; i++) {
        for (j = rect->LefttopPoint.col; j < rect->LefttopPoint.col + rect->width; j++) {
            if ((int) GradientImg->DirectionImage->data[i][j] % AngularResolution > AngularResolution / 2) {
                temp = ((int) GradientImg->DirectionImage->data[i][j] / AngularResolution + 1) * AngularResolution;
                index = temp / AngularResolution;
                Gradhist->norm[index] = Gradhist->norm[index] + GradientImg->NormImage->data[i][j];
            } else {
                temp = (int) GradientImg->DirectionImage->data[i][j] / AngularResolution * AngularResolution;
                index = temp / AngularResolution;
                Gradhist->norm[index] = Gradhist->norm[index] + GradientImg->NormImage->data[i][j];
            }
        }
    }

    /*
    AngleNumber = 360 / AngularResolution;
    height = Img->height;
    width = Img->width;

    Gradhist = (GradHist*)malloc(sizeof(GradHist));
    Gradhist->angle = (int*)malloc(sizeof(int)* AngleNumber);
    Gradhist->norm = (double*)malloc(sizeof(double)* AngleNumber);

    for (i = 0; i < AngleNumber; i++){
    Gradhist->angle[i] = i * AngularResolution;
    }

    for (i = 1; i < width - 1; i++){
    for (j = 1; j < height - 1; j++){
    point.column = i;
    point.row = j;
    grad = GetGradient(&point, Img);
    if ((int)grad->direction % AngularResolution > AngularResolution / 2){
    grad->direction = ((int)grad->direction / AngularResolution + 1) * AngularResolution;
    index = grad->direction / AngularResolution - 1;
    Gradhist->norm[index] = grad->norm;
    }
    else
    {
    grad->direction = (int)grad->direction / AngularResolution * AngularResolution;
    index = grad->direction / AngularResolution - 1;
    Gradhist->norm[index] = grad->norm;
    }
    }
    }
    */

    return Gradhist;
}


double *Normalize(double *Array, int length) {
    double *NormalizedArray;
    double sum;
    int i;

    NormalizedArray = new double[length];

    sum = 0.0;

    for (i = 0; i < length; i++) {
        sum = sum + Array[i];
    }

    for (i = 0; i < length; i++) {
        NormalizedArray[i] = Array[i] / sum;
    }

    return NormalizedArray;
}


void CreateKeyPointDescriptor(ImagePyramid *DiffPyramid, ImagePyramid *ImgPyramid, KeyPointGroup *KeyPG) {
    GradHist *TempGradHist;
    Rect Rect;
    Matrix *RegionImg;
    Matrix *RotateImg;
    GradientImage *GradientImg;
    double *Tempdesc;
    int KeyPointNumber;
    int RegionSize;
    int SquareSize;
    int MainDerection;
    int AngularResolution;
    int width;
    int height;
    int radius;
    int i;
    int j;
    int m;

    KeyPointNumber = KeyPG->number;
    AngularResolution = 45;

    //启用所有特征点
    for (i = 0; i < KeyPointNumber; i++) {
        KeyPG->point[i].flag = true;
    }

    for (i = 1; i < KeyPointNumber; i++) {
        RegionSize = 5 * 3 * KeyPG->point[i].sigma;
        SquareSize = RegionSize / 5;
        radius = 0.71 * RegionSize;

        //如果特征点的邻域超过边界，则设置flag为false，不计算描述子
        width = ImgPyramid->mats[KeyPG->point[i].octavenum][KeyPG->point[i].intervalnum].width;
        height = ImgPyramid->mats[KeyPG->point[i].octavenum][KeyPG->point[i].intervalnum].height;
        if (KeyPG->point[i].point.row - radius < 0 || KeyPG->point[i].point.col - radius < 0
            || KeyPG->point[i].point.col + radius > width || KeyPG->point[i].point.row + radius > height) {
            KeyPG->point[i].flag = false;
            continue;
        }

        KeyPG->point[i].descriptor = (double *) malloc(sizeof(double) * 128);
        MainDerection = KeyPG->point[i].grad.direction;
        Rect.LefttopPoint.row = KeyPG->point[i].point.row - radius;//特征向量一行一行的排序
        Rect.LefttopPoint.col = KeyPG->point[i].point.col - radius;
        Rect.height = 2 * radius;
        Rect.width = 2 * radius;

        RegionImg = CropImage(&(ImgPyramid->mats[KeyPG->point[i].octavenum][KeyPG->point[i].intervalnum]), &Rect);
        RotateImg = RotateImage(RegionImg, MainDerection);
        GradientImg = GetGradientImage(RotateImg);

        //debug用
        /*
        ClImage* img_result1;
        img_result1 = SaveImage(RegionImg);
        bool flag1 = clSaveImage("C:/Users/Administrator/Desktop/region.bmp", img_result1);
        img_result1 = SaveImage(RotateImg);
        flag1 = clSaveImage("C:/Users/Administrator/Desktop/ratate.bmp", img_result1);
        */

        for (j = 0; j < 16; j++) {
            Rect.LefttopPoint.row = radius - SquareSize * 2 + (j % 4) * SquareSize;//特征向量一行一行的排序
            Rect.LefttopPoint.col = radius - SquareSize * 2 + (j / 4) * SquareSize;
            Rect.height = SquareSize;
            Rect.width = SquareSize;

            TempGradHist = GetGradHist(GradientImg, AngularResolution, &Rect);

            for (m = 0; m < 8; m++) {
                KeyPG->point[i].descriptor[j * 8 + m] = TempGradHist->norm[m];
            }

        }

        //归一化
        Tempdesc = Normalize(KeyPG->point[i].descriptor, 128);
        for (j = 0; j < 128; j++) {
            KeyPG->point[i].descriptor[j] = Tempdesc[j];
        }
    }
}


KeyPointGroup *SIFT(Matrix *mat) {

    ImagePyramid *ip = new ImagePyramid;
    std::cout << "build GuassPyramif" << std::endl;
    ip->build(mat, 64, 1);
    std::cout << "build GuassDiffPyramif" << std::endl;
    ImagePyramid *diff_ip = ip->buildDiff();
    std::cout << "SearchKeyPointPosition" << std::endl;

    auto KP = SearchKeyPointPosition(ip, diff_ip);
    std::cout << "CreateKeyPointDescriptor" << std::endl;
    CreateKeyPointDescriptor(diff_ip, ip, KP);
    return KP;
}


Point GetInitPosition(KeyPoint *KP) {
    Point p;

    p.row = KP->point.row * pow((double) 2, KP->octavenum);
    p.col = KP->point.col * pow((double) 2, KP->octavenum);

    return p;
}


double Norm(double *p, double *q, int length) {
    double Answer;
    double sum;
    int i;

    sum = 0;
    for (i = 0; i < length; i++) {
        sum = sum + pow(p[i] - q[i], 2);
    }

    Answer = pow(sum, 0.5);

    return Answer;
}


MatchPoint *Match(KeyPointGroup *KPG1, KeyPointGroup *KPG2) {
    MatchPoint *Matchpoint;
    int KPG1Number;
    int KPG2Number;
    int MatchNumber;
    int i, j;
    double nn1, nn2;
    double TempNorm;
    int index;
    double threshold;

    Matchpoint = (MatchPoint *) malloc(sizeof(MatchPoint));
    Matchpoint->Img1 = (Point *) malloc(sizeof(Point));
    Matchpoint->Img2 = (Point *) malloc(sizeof(Point));
    Matchpoint->Matchnumber = 0;
    nn1 = DBL_MAX;
    nn2 = DBL_MAX;
    threshold = 0.9;
    MatchNumber = 0;

    KPG1Number = KPG1->number;
    KPG2Number = KPG2->number;

    for (i = 1; i < KPG1Number; i++) {
        if (KPG1->point[i].flag == false) {
            continue;
        }

        for (j = 1; j < KPG2Number; j++) {
            if (KPG2->point[j].flag == false) {
                continue;
            }
            /*
             * Norm是两个特征向量之间的差向量的模，模越小，说明两个向量越相似，以此来匹配特征点
             */
            TempNorm = Norm(KPG1->point[i].descriptor, KPG2->point[j].descriptor, 128);
            if (TempNorm < nn1) {
                nn1 = TempNorm;
                index = j;
            } else if (TempNorm < nn2) {
                nn2 = TempNorm;
            }
        }

        /*
         * nn1表示相似性的最小值，nn2表示相似性的第二小值，
         * 如果nn1与nn1之间很接近的话，nn1和nn2对应的两个待匹配点就很难区分，这两个点会互相干扰匹配，
         * 所以要判断nn1和nn2之间的相似性来最终确定匹配点
         */
        if (nn1 / nn2 < threshold) {
            MatchNumber++;
            Matchpoint->Img1 = (Point *) realloc(Matchpoint->Img1, sizeof(Point) * MatchNumber);
            Matchpoint->Img2 = (Point *) realloc(Matchpoint->Img2, sizeof(Point) * MatchNumber);
            Matchpoint->Img1[MatchNumber - 1] = GetInitPosition(&KPG1->point[i]);
            Matchpoint->Img2[MatchNumber - 1] = GetInitPosition(&KPG2->point[index]);
            //如果某个特征点匹配，则不再参与匹配
            KPG2->point[index].flag = false;
            nn1 = DBL_MAX;
            nn2 = DBL_MAX;
        }
    }

    Matchpoint->Matchnumber = MatchNumber;

    return Matchpoint;
}


Matrix *DrawCircle(Matrix *Img, Point *p, int size) {
    if (p->col + (size - 1) / 2 >= Img->width || p->col - (size - 1) / 2 < 0
        || p->row + (size - 1) / 2 >= Img->height || p->row - (size - 1) / 2 < 0) {
        return Img;
    }

    double r;
    double tempr;
    int i, j;

    r = size / 2;

    for (i = -(size - 1) / 2; i <= (size - 1) / 2; i++) {
        for (j = -(size - 1) / 2; j <= (size - 1) / 2; j++) {
            tempr = pow((double) i, 2) + pow((double) j, 2);
            tempr = pow(tempr, 0.5);
            if (fabs((tempr - r)) < 0.5) {
                Img->data[p->row + i][p->col + j] = 255;
            }
        }
    }
}


Matrix *CreateMatchMatrix(Matrix *Img1, Matrix *Img2, MatchPoint *Matchpoint) {
    Matrix *MatchImg;
    int linewidth = 2;
    Point *tempp;
    int height;
    int width;
    int i, j, m;
    double k;
    double x1, y1, x2, y2;

    tempp = new Point;
    MatchImg = new Matrix;

    height = Img1->height;
    width = Img1->width + Img2->width;
    MatchImg->height = height;
    MatchImg->width = width;

    MatchImg->data = new double *[height];
    for (i = 0; i < height; i++) {
        MatchImg->data[i] = new double[width];
    }

    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
            if (j < Img1->width) {
                MatchImg->data[i][j] = Img1->data[i][j];
                MatchImg->data[i][j] = Img1->data[i][j];
                MatchImg->data[i][j] = Img1->data[i][j];
            } else {
                MatchImg->data[i][j] = Img2->data[i][j - Img1->width];
                MatchImg->data[i][j] = Img2->data[i][j - Img1->width];
                MatchImg->data[i][j] = Img2->data[i][j - Img1->width];
            }
        }
    }

    for (m = 0; m < Matchpoint->Matchnumber; m++) {
        x1 = Matchpoint->Img1[m].col;
        y1 = Matchpoint->Img1[m].row;
        x2 = Img1->width + Matchpoint->Img2[m].col;
        y2 = Matchpoint->Img2[m].row;
        k = (y2 - y1) / (x2 - x1);

        for (j = 0; j < Matchpoint->Img2[m].col + Img1->width - Matchpoint->Img1[m].col; j++) {
            i = Matchpoint->Img1[m].row + k * j;
            tempp->col = Matchpoint->Img2[m].col + Img1->width;
            tempp->row = Matchpoint->Img2[m].row;

            DrawCircle(MatchImg, &Matchpoint->Img1[m], 7);
            DrawCircle(MatchImg, tempp, 7);

            MatchImg->data[i][Matchpoint->Img1[m].col + j] = 255;
            MatchImg->data[i + 1][Matchpoint->Img1[m].col + j] = 255;
        }
    }

    return MatchImg;
}

#endif //_DIFT_HPP_