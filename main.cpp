#include <iostream>
#include <fstream>
#include <cstring>
#include <windows.h>
#include <zconf.h>
#include "utils.h"
#include "sift.hpp"
#include "pgmer.hpp"
#include "img.h"

void DrawImg(Matrix *mat, int x, int y){
    HWND wnd;	//窗口句柄
    HDC dc;	//绘图设备环境句柄
    wnd = GetForegroundWindow(); //获取窗口句柄
    dc = GetDC(wnd);	//获取绘图设备
    for(int j=0;j< mat->height;j++) {
        for (int i = 0; i < mat->width; i++) {
            float v = mat->data[j][i];
            SetPixel(dc, i + x, j + y, RGB(v,v,v)); //画像素点
        }
    }
    std::cout << "any key next..." << std::endl;
    getchar();
}


int main(){
    Pgmer p1;
    p1.ReadImg("01.pgm");
    std::cout << "convert img 1 to mat" << std::endl;
    Matrix *mat1 = p1.To2DMatrix();
    KeyPointGroup *kp1 =  SIFT(mat1);

    std::cout << "img 1 have " << kp1->number - 1 << " keypoint" << std::endl;

    Pgmer p2;
    p2.ReadImg("02.pgm");
    std::cout << "convert img 2 to mat" << std::endl;
    Matrix *mat2 = p2.To2DMatrix();
    KeyPointGroup *kp2 =  SIFT(mat2);
    std::cout << "img 2 have " << kp2->number - 1 << " keypoint" << std::endl;


    std::cout << "matching tow img..." << std::endl;
    MatchPoint *matchPoint = Match(kp1, kp2);
    std::cout << "match point:" << matchPoint->Matchnumber << std::endl;

    std::cout << "creatwMatchImg..." << std::endl;
    auto result = CreateMatchMatrix(mat1, mat2, matchPoint);

    std::cout << "generate img..." << std::endl;
    result->WriteImg("result.pgm");
    DrawImg(result, 10, 10);
    std::cout << "done" << std::endl;

    return 0;
}