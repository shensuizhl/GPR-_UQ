1.注意改变：输入观测数据文件及观测值数量，似然函数，模拟值的计算与输入
how to run makefile（gfortran中需要）
2.dreamtest.txt文件，注意：napai<nchain
3.input_point.txt 观测点时空位置文件//格式 x y t 其中第一行是in_dim(表示坐标的维数如（x,y,t）为3)
4.predict_point.txt预测点时空位置文件//格式 x y t 其中第一行是dim_pre ,seed,api
对于输入格式200，190要进行修改
5.prior.in文件，注意测量误差一项，尽量不为0,否则会出现非对称正定矩阵

