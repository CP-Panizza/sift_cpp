#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_


#include <fstream>

class Matrix {
public:
    int width;
    int height;
    double **data = nullptr;

    Matrix() {}

    ~Matrix(){
        for (int i = 0; i < this->height; ++i) {
            delete[](this->data[i]);
        }
        delete[](this->data);
    }

    Matrix(int h, int w, int *_data) {
        this->data = new double *[h];
        for (int i = 0; i < h; ++i) {
            this->data[i] = new double[w];
        }

        for (int j = 0; j < h; ++j) {
            for (int i = 0; i < w; ++i) {
                this->data[j][i] = _data[w * j + i];
            }
        }

        this->width = w;
        this->height = h;
    }

    Matrix *Copy(){
        auto t = new Matrix;
        t->width = this->width;
        t->height = this->height;
        t->data = new double *[this->height];
        for (int i = 0; i < this->height; ++i) {
            t->data[i] = new double[this->width];
        }

        for (int j = 0; j < this->height; ++j) {
            for (int i = 0; i < this->width; ++i) {
                t->SetPix(j, i, this->GetPix(j, i));
            }
        }
        return t;
    }

    void WriteImg(const std::string name) {
        std::ofstream f(name);
        if (f.is_open()) {
            f << "P2\n" << this->width << " " << this->height << "\n255\n";
            int k = 1;
            for (int i = 0; i < this->height; ++i) {
                for (int j = 0; j < this->width; ++j) {
                    if (k % 10) {
                        f << this->data[i][j] << " ";
                    } else {
                        f << this->data[i][j] << "\n";
                    }
                    k++;
                }
            }
        } else {
            perror("open outfile err!");
        }
        f.close();
    }

    bool SetPix(int row, int col, float val) {
        if (col > this->width || row > this->height) return false;
        this->data[row][col] = val;
        return true;
    }

    float GetPix(int row, int col) {
        if (col > this->width || row > this->height) return -1;
        return this->data[row][col];
    }

    Matrix *DownSampling() {
        Matrix *down = new Matrix;

        down->height = this->height / 2;
        down->width = this->width / 2;

        down->data = new double*[down->height];
        for (int i = 0; i < down->height; ++i) {
            down->data[i] = new double[down->width];
        }

        for (int i = 0; i < down->height; ++i) {
            for (int j = 0; j < down->width; ++j) {
                down->SetPix(i, j, this->GetPix(i * 2, j * 2));
            }
        }

        return down;
    }

};


#endif