//
//  NISwGSP_Stitching.cpp
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#include "NISwGSP_Stitching.h"

NISwGSP_Stitching::NISwGSP_Stitching(const MultiImages & _multi_images) : MeshOptimization(_multi_images) {
    
}

void NISwGSP_Stitching::setWeightToAlignmentTerm(const double _weight) {
    MeshOptimization::setWeightToAlignmentTerm(_weight);
}

void NISwGSP_Stitching::setWeightToLocalSimilarityTerm(const double _weight) {
    MeshOptimization::setWeightToLocalSimilarityTerm(_weight);
}

void NISwGSP_Stitching::setWeightToGlobalSimilarityTerm(const double _weight_beta,
                                                        const double _weight_gamma,
                                                        const enum GLOBAL_ROTATION_METHODS _global_rotation_method) {
    MeshOptimization::setWeightToGlobalSimilarityTerm(_weight_beta, _weight_gamma, _global_rotation_method);
}

void NISwGSP_Stitching::setWeightToLinePreserveTerm(const double _weight) {
    MeshOptimization::setWeightToLinePreserveTerm(_weight);
}

// solve应分为两阶段求解，第一阶段求解目标为对齐、保相似性，第二阶段求解目标在第一阶段基础上增加保直线
Mat NISwGSP_Stitching::solve(const BLENDING_METHODS & _blend_method) {
    const MultiImages & multi_images = getMultiImages();
    
    vector<Triplet<double> > triplets;      // triplets（三元数）代表优化项，存储优化项在矩阵中的行列坐标和值，即稀疏矩阵的各项
    vector<pair<int, double> > b_vector;    // b_vector（pair<int, double>)代表Ax=b右向量，存储行坐标和值即优化目标，
    
    // 第一阶段优化
    cout << "*************第一阶段优化**************" << endl;
    const pair<int, int> & origin_sparse_size = reserveData(triplets, b_vector, DIMENSION_2D);
    
    triplets.emplace_back(0, 0, STRONG_CONSTRAINT);
    triplets.emplace_back(1, 1, STRONG_CONSTRAINT);
    b_vector.emplace_back(0,    STRONG_CONSTRAINT);
    b_vector.emplace_back(1,    STRONG_CONSTRAINT);
    
    prepareAlignmentTerm(triplets); // 实际上是用feature pairs而不是用的mesh pairs
    prepareSimilarityTerm(triplets, b_vector);  // 准备全局相似项和局部相似项
    //prepareLinePreserveTerm(triplets, b_vector); // 准备全局直线结构保护箱
    
    vector<vector<Point2> > new_vertices;
    
    new_vertices = getImageVerticesBySolving(triplets, b_vector, true);  // 解稀疏矩阵，求最优点集
    
#ifndef NDEBUG
    // 输出优化后的顶点集
    ofstream outFile;
    outFile.open(multi_images.parameter.temp_dir + multi_images.parameter.file_name + "-line_control_points.txt");
//    for (int i = 0; i < new_vertices.size(); i++) {
//        outFile << "第" << i << "张图片计算顶点：" << endl;
//        for (int j = 0; j < multi_images.images_data[i].mesh_2d->nh; j++) {
//            for (int k = 0; k < multi_images.images_data[i].mesh_2d->nw; k++) {
//                outFile << new_vertices[i][j*multi_images.images_data[i].mesh_2d->nw+k] << " ";
//            }
//            outFile << endl;
//        }
//    }
    // 打印直线控制点端点位置
    outFile << "****************第一阶段优化***************" << endl;
    const vector<vector<vector<InterpolateVertex> > > & interpolateVerticesOfSelectedLines = multi_images.getInterpolateVerticesOfSelectedLines();
    for (int i = 0; i < multi_images.images_data.size(); i++) {
        outFile << "第" << i << "张图片：" << endl;
        for (int j = 0; j < interpolateVerticesOfSelectedLines[i].size(); j++) {
            outFile << "第" << j << "条直线控制端点位置: " << endl;
            for (int p = 0; p < interpolateVerticesOfSelectedLines[i][j].size(); p++) {
                outFile << multi_images.images_data[i].mesh_2d->getPointFromInterpolateVertex(interpolateVerticesOfSelectedLines[i][j][p], new_vertices[i]) << " ";
            }
            outFile << endl;
        }
    }

    outFile.close();
#endif
    
    // 第二阶段优化
    cout << "*************第二阶段优化**************" << endl;
    reserveExtraData(triplets, b_vector, origin_sparse_size);
    
    prepareLinePreserveTerm(triplets, b_vector, new_vertices);    // 准备全局直线结构保护项
    
    vector<vector<Point2> > final_vertices;
    
    final_vertices = getImageVerticesBySolving(triplets, b_vector, false);
    
#ifndef NDEBUG
    outFile.open(multi_images.parameter.temp_dir + multi_images.parameter.file_name + "-line_control_points.txt", ios::app);
    // 打印直线控制点端点位置
    outFile << "****************第二阶段优化****************" << endl;
    for (int i = 0; i < multi_images.images_data.size(); i++) {
        outFile << "第" << i << "张图片：" << endl;
        for (int j = 0; j < interpolateVerticesOfSelectedLines[i].size(); j++) {
            outFile << "第" << j << "条直线控制端点位置: " << endl;
            for (int p = 0; p < interpolateVerticesOfSelectedLines[i][j].size(); p++) {
                outFile << multi_images.images_data[i].mesh_2d->getPointFromInterpolateVertex(interpolateVerticesOfSelectedLines[i][j][p], final_vertices[i]) << " ";
            }
            outFile << endl;
        }
    }
    
    outFile.close();
#endif
    
    Size2 target_size = normalizeVertices(final_vertices);   // 归一化点集
    
    Mat result = multi_images.textureMapping(final_vertices, target_size, _blend_method);
#ifndef NDEBUG
    multi_images.writeResultWithMesh(result, final_vertices, "-[NISwGSP]" +
                                     GLOBAL_ROTATION_METHODS_NAME[getGlobalRotationMethod()] +
                                     BLENDING_METHODS_NAME[_blend_method] +
                                     "[Mesh]", false);
    multi_images.writeResultWithMesh(result, final_vertices, "-[NISwGSP]" +
                                     GLOBAL_ROTATION_METHODS_NAME[getGlobalRotationMethod()] +
                                     BLENDING_METHODS_NAME[_blend_method] +
                                     "[Border]", true);
    if (getLinePreserveTermWeight()) {
        multi_images.writeResultWithLines(result, final_vertices, "-[NISwGSP]" +
                                          GLOBAL_ROTATION_METHODS_NAME[getGlobalRotationMethod()] +
                                          BLENDING_METHODS_NAME[_blend_method] +
                                          "[LINES]");
    }
#endif
    return result;
}

void NISwGSP_Stitching::writeImage(const Mat & _image, const string _post_name) const {
    const MultiImages & multi_images = getMultiImages();
    const Parameter & parameter = multi_images.parameter;
    string file_name = parameter.file_name;
    
    imwrite(parameter.result_dir + file_name + "-" +
            "[NISwGSP]" +
            GLOBAL_ROTATION_METHODS_NAME[getGlobalRotationMethod()] +
            _post_name +
            ".png", _image);
}