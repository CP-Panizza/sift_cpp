#ifndef _GAUSSFILTER_HPP_
#define _GAUSSFILTER_HPP_

//高斯卷积核
class GaussFilter {
public:
    int Size;
    float **FilterData;

    GaussFilter(double Sigma){
        double sum = 0.0;
        auto Size = static_cast<int>(6 * Sigma + 1);
        if(Size % 2 == 0){
            Size += 1;
        }
        this->Size = Size;

        this->FilterData = new float *[this->Size];
        for (int i = 0; i < this->Size; ++i) {
            this->FilterData[i] = new float[this->Size];
        }

        for (int i = 0; i < Size; i++){
            for (int j = 0; j < Size; j++){
                this->FilterData[i][j] = static_cast<float>(exp(-((i - Size / 2) * (i - Size / 2) + (j - Size / 2) * (j - Size / 2)) / (2.0 * Sigma * Sigma)));
                sum += this->FilterData[i][j];
            }
        }

        for (int i = 0; i < Size; i++){
            for (int j = 0; j < Size; j++){
                this->FilterData[i][j] = static_cast<float>(this->FilterData[i][j] / sum);
            }
        }
    }

    ~GaussFilter(){
        for (int i = 0; i < this->Size; ++i) {
            delete[](this->FilterData[i]);
        }
        delete[](this->FilterData);
    }

    /**
     * 对矩阵进行卷积，返回卷积后的矩阵
     * @param data
     * @return
     */
    Matrix *Convolution(Matrix *data){
        float sum = 0;
        int bound = (this->Size - 1) / 2;
        Matrix *conv = new Matrix;
        conv->width = data->width;
        conv->height = data->height;
        conv->data = new double*[conv->height];
        for (int i = 0; i < conv->height; ++i) {
            conv->data[i] = new double[conv->width];
        }
        Matrix *padding_data = this->Padding(data);

        for (int j = 0; j < data->height; ++j) {
            for (int i = 0; i < data->width; ++i) {
                for (int k = 0; k < this->Size; ++k) {
                    for (int m = 0; m < this->Size; ++m) {
                        sum += padding_data->data[j + k][i + m] * this->FilterData[k][m];
                    }
                }
                conv->data[j][i] = sum; //转成int可方便导出pgm查看
                sum = 0;
            }
        }
        delete(padding_data);
        return conv;
    }



    /**
     * 在卷积之前对矩阵预先进行边界padding
     * @param data
     * @return
     */

    Matrix *Padding(Matrix *data){
        Matrix *padding = new Matrix;
        int bound = (this->Size - 1) / 2;
        padding->width = data->width + this->Size - 1;
        padding->height = data->height + this->Size - 1;

        padding->data = new double *[padding->height];
        for (int i = 0; i < padding->height; ++i) {
            padding->data[i] = new double[padding->width];
        }

        for (int j = 0; j < padding->height; ++j) {
            for (int i = 0; i < padding->width; ++i) {
                if(j > (bound - 1) && j <= (padding->height - bound - 1) && i > (bound - 1) && i <=(padding->width - bound - 1)){
                    padding->data[j][i] = data->data[j - bound][i - bound];
                } else {
                    padding->data[j][i] = 0;
                }
            }
        }
        return padding;
    }
};




#endif //_GAUSSFILTER_HPP_