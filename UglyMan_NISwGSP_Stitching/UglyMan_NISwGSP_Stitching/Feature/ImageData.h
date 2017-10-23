//
//  ImageData.h
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#ifndef __UglyMan_Stitiching__ImageData__
#define __UglyMan_Stitiching__ImageData__

#include <memory>
#include "Statistics.h"
#include "StringTools.h"
#include "FeatureController.h"
#include "MeshGrid.h"

void onMouse(int event, int x, int y, int flags, void *param);

struct UserData {
    bool isDelete = false;
    RNG g_rng;
    Point pre, cur;
    Mat img, tmp;
    vector<Vec4d> lines;
};

class LineData {
public:
    LineData(const Point2 & _a,
             const Point2 & _b,
             const double _width,
             const double _length);
    Point2 data[2];
    double width, length;
private:
};

class LineSegments {
public:
    LineSegments(const vector<Point2> & _points,
                 const double _length,
                 const double _step);
    vector<Point2> points;
    double length, step;
private:
};

// 这里非常巧妙的应用了typedef替代重复的声明过程
typedef const bool (LINES_FILTER_FUNC)(const double _data, \
                                       const Statistics & _statistics);

LINES_FILTER_FUNC LINES_FILTER_NONE;
LINES_FILTER_FUNC LINES_FILTER_WIDTH;
LINES_FILTER_FUNC LINES_FILTER_LENGTH;


class ImageData {
public:
    string file_name, file_extension;
    const string * file_dir, * debug_dir, * temp_dir;
    ImageData(const string & _file_dir,
              const string & _file_full_name,
              LINES_FILTER_FUNC * _width_filter,
              LINES_FILTER_FUNC * _length_filter,
              const string * _debug_dir = NULL,
              const string * _temp_dir = NULL);
    
    const Mat & getGreyImage() const;
    const vector<LineData> & getLines() const;
    // 选择需要的直线原始点集
    const vector<LineSegments> & getSelectedLines() const;
    const vector<Point2> & getFeaturePoints() const;
    const vector<FeatureDescriptor> & getFeatureDescriptors() const;
    
    void clear();
    
    Mat img, rgba_img, alpha_mask;
    unique_ptr<Mesh2D> mesh_2d;
    
private:
    LINES_FILTER_FUNC * width_filter, * length_filter;
    
    mutable Mat grey_img;
    mutable vector<LineData> img_lines;
    // 选中直线的点集
    mutable vector<LineSegments> selected_lines;
    mutable vector<Point2> feature_points;
    mutable vector<FeatureDescriptor> feature_descriptors;
};

#endif /* defined(__UglyMan_Stitiching__ImageData__) */
