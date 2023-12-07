English Version:

Note that you need to download MATLAB Runtime R2018b in advance.

The types of files that need to be placed in the folder：
1. Point cloud rotated according to the direction of the sun, with the suffix .log and the format x, y, z
2. Original point cloud with normal vector, suffix is .​​txt, format is x, y, z, nx, ny, nz
3. The three-dimensional distribution of the scattering penetration coefficient. The file is named GapFractionVertical_Mean.txt. It is a matrix of N*1, and N is the number of scattering layers =fix[(max_height-min_height)/height_intervel]+1;
4. Incident solar radiation and direction. The file name is incident_sun.txt, which is a 1*5 matrix. The five columns are: solar incident direct PAR, solar incident scattered PAR, and solar incident vectors s_nx, s_ny, s_nz.
5. The parameters of the model are a 1*3 matrix, including the large voxel size and small voxel size when calculating direct radiation, and the height interval for calculating scattering.

How to operate:
1. Double-click rad_dir_dif.exe and select the point cloud with the suffix .log
2. Select the point cloud with the normal vector suffix txt
3. Just wait. The output file format is rad*.log, and the format is x, y, z, direct PAR, scattered PAR, and total PAR;


中文版：
注意需要提前下载 MATLAB Runtime R2018b。 

文件夹里面需要放置的文件类型

1、根据太阳方向旋转过的点云，后缀为.log，格式为x，y，z
2、有法向量的原始数据点云，后缀为txt，格式为x，y，z，nx，ny，nz
3、散射穿透系数的三维分布，文件命名为GapFractionVertical_Mean.txt，是N*1的矩阵，N为散射分层数=fix[（max_height-min_height）/height_intervel]+1;
4、入射太阳辐射以及方向，文件名命名为incident_sun.txt，是1*5的矩阵，5列分别为：太阳入射直射PAR、太阳入射散射PAR，以及太阳入射向量s_nx，s_ny，s_nz。
5、模型的参数，是1*3的矩阵，包括计算直射时的大体元大小、小体元大小，以及计算散射的高度间隔

操作方法：
1、双击rad_dir_dif.exe，选取后缀为.log的点云
2、选取有法向量后缀为txt的点云
3、等待即可。输出文件格式为rad*.log，格式为x，y，z，直射PAR、散射PAR、总PAR；


