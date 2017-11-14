//
//  ImageData.cpp
//  UglyMan_Stitching
//
//  Created by uglyman.nothinglo on 2015/8/15.
//  Copyright (c) 2015 nothinglo. All rights reserved.
//

#include "ImageData.h"

LineData::LineData(const Point2 & _a,
                   const Point2 & _b,
                   const double _width,
                   const double _length) {
    data[0] = _a;
    data[1] = _b;
    width   = _width;
    length  = _length;
}

const bool LINES_FILTER_NONE(const double _data,
                             const Statistics & _statistics) {
    return true;
};

const bool LINES_FILTER_WIDTH (const double _data,
                               const Statistics & _statistics) {
    return _data >= MAX(2.f, (_statistics.min + _statistics.mean) / 2.f);
    return true;
};

const bool LINES_FILTER_LENGTH(const double _data,
                               const Statistics & _statistics) {
    return _data >= MAX(10.f, _statistics.mean);
    return true;
};


ImageData::ImageData(const string & _file_dir,
                     const string & _file_full_name,
                     LINES_FILTER_FUNC * _width_filter,
                     LINES_FILTER_FUNC * _length_filter,
                     const string * _debug_dir,
                     const string * _temp_dir) {
    
    file_dir = &_file_dir;
    std::size_t found = _file_full_name.find_last_of(".");
    assert(found != std::string::npos);
    file_name = _file_full_name.substr(0, found);
    file_extension = _file_full_name.substr(found);
    debug_dir = _debug_dir;
    temp_dir = _temp_dir;

    grey_img = Mat();
    
    width_filter  = _width_filter;
    length_filter = _length_filter;
    
    img = imread(*file_dir + file_name + file_extension);
    rgba_img = imread(*file_dir + file_name + file_extension, IMREAD_UNCHANGED);
    
    float original_img_size = img.rows * img.cols;
    
    if(original_img_size > DOWN_SAMPLE_IMAGE_SIZE) {
        float scale = sqrt(DOWN_SAMPLE_IMAGE_SIZE / original_img_size);
        resize(img, img, Size(), scale, scale);
        resize(rgba_img, rgba_img, Size(), scale, scale);
    }
    
    assert(rgba_img.channels() >= 3);
    if(rgba_img.channels() == 3) {
        cvtColor(rgba_img, rgba_img, CV_BGR2BGRA);
    }
    vector<Mat> channels;
    split(rgba_img, channels);
    alpha_mask = channels[3];
    mesh_2d = make_unique<MeshGrid>(img.cols, img.rows);
}

// 转为灰度图矩阵
const Mat & ImageData::getGreyImage() const {
    if(grey_img.empty()) {
        cvtColor(img, grey_img, CV_BGR2GRAY);
    }
    return grey_img;
}

// 获取图片的LSD特征，也是用来作为直线结构保持项的
const vector<LineData> & ImageData::getLines() const {
    if(img_lines.empty()) {
        const Mat & grey_image = getGreyImage();
        Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_STD);
        
        vector<Vec4f>  lines;
        vector<double> lines_width, lines_prec, lines_nfa;
        ls->detect(grey_image, lines, lines_width, lines_prec, lines_nfa);
        
        vector<double> lines_length;
        vector<Point2> lines_points[2];
        
        const int line_count = (int)lines.size();
        
        lines_length.reserve(line_count);
        lines_points[0].reserve(line_count);
        lines_points[1].reserve(line_count);
        
        for(int i = 0; i < line_count; ++i) {
            lines_points[0].emplace_back(lines[i][0], lines[i][1]);
            lines_points[1].emplace_back(lines[i][2], lines[i][3]);
            lines_length.emplace_back(norm(lines_points[1][i] - lines_points[0][i]));
        }
        
        const Statistics width_statistics(lines_width), length_statistics(lines_length);
        for(int i = 0; i < line_count; ++i) {
            if( width_filter( lines_width[i],  width_statistics) &&
               length_filter(lines_length[i], length_statistics)) {
                img_lines.emplace_back(lines_points[0][i],
                                       lines_points[1][i],
                                       lines_width[i],
                                       lines_length[i]);
            }
        }
#ifndef NDEBUG
        vector<Vec4f> draw_lines;
        draw_lines.reserve(img_lines.size());
        for(int i = 0; i < img_lines.size(); ++i) {
            draw_lines.emplace_back(img_lines[i].data[0].x, img_lines[i].data[0].y,
                                    img_lines[i].data[1].x, img_lines[i].data[1].y);
        }
        Mat canvas = Mat::zeros(grey_image.rows, grey_image.cols, grey_image.type());
        ls->drawSegments(canvas, draw_lines);
        imwrite(*debug_dir + "line-result-" + file_name + file_extension, canvas);
#endif
    }
    return img_lines;
}

// ***得到选定的直线 自动化处理方式暂放置一边
//const vector<LineSegments> & ImageData::getSelectedLines() const {
//    if(selected_lines.empty()) {
//        const vector<LineData> & lines = getLines();
//        for(int i = 0; i < lines.size(); i++) {
//            // 选择较长的直线，在该直线上采样
//            if(lines[i].length >= LINES_THRESHOLD) {
//                vector<Point2> line_segments;
//                double step = 10.0;
//                double stepX = step * (lines[i].data[1].x - lines[i].data[0].x)/lines[i].length;
//                double stepY = step * (lines[i].data[1].y - lines[i].data[0].y)/lines[i].length;
//                for(double x = lines[i].data[0].x, y = lines[i].data[0].y; (x - lines[i].data[1].x) * stepX < 0; x += stepX, y += stepY) {
//                    line_segments.emplace_back(x, y);
//                }
//                line_segments.emplace_back(lines[i].data[1]);
//                selected_lines.emplace_back(line_segments, lines[i].length, step);
//            }
//        }
//#ifndef NDEBUG
//        vector<Vec4f> draw_lines;
//        const Mat & grey_image = getGreyImage();
//        Mat canvas = Mat::zeros(grey_image.rows, grey_image.cols, grey_image.type());
//        draw_lines.reserve(selected_lines.size());
//        for(int i = 0; i < selected_lines.size(); i++) {
//            line(canvas, selected_lines[i].points[0], selected_lines[i].points[selected_lines[i].points.size()-1],
//                 Scalar(0, 0, 255), 1, LINE_AA);
//        }
//        imwrite(*debug_dir + "selected-lines-" + file_name + file_extension, canvas);
//#endif
//    }
//    return selected_lines;
//}

// 得到长直线结果集 TODO : 寻求自动化提取长直线结构的方式
const vector<vector<Point2> > & ImageData::getSelectedLines() const {
    if(selected_lines.empty()) {
        // choose to load the stored lines or select by gui
        // if load the lines, don't need to store them again
        // else, select by gui and store them
        vector<Vec4d> lines;

        cout << "load the lines data? (y/n)" << endl;
        if (cin.get() == 'y') {
            cin.get();
            ifstream inFile(*temp_dir + file_name + "_selectedlines.txt");
            if(!inFile.is_open()) {
                cout << "Error opening file" << endl;
                exit(0);
            } else {
                string str;
                while (getline(inFile, str)) {
                    const vector<string> & strs = split(str.substr(1, str.size()-2), ", ");
                    lines.emplace_back(atof(strs[0].c_str()),
                                       atof(strs[1].c_str()),
                                       atof(strs[2].c_str()),
                                       atof(strs[3].c_str()));
                }
                inFile.close();
            }
        } else {
            UserData userdata;
            img.copyTo(userdata.img);
            img.copyTo(userdata.tmp);
            
            namedWindow("selectLines");
            setMouseCallback("selectLines", onMouse, (void*)&userdata);
            
            while (1) {
                imshow("selectLines", userdata.img);
                if(waitKey(10) == 27) {
                    break;
                }
            }
            destroyAllWindows();
            
            imwrite(*debug_dir + "selected-lines-" + file_name + file_extension, userdata.img);
            
            ofstream outFile(*temp_dir + file_name + "_selectedlines.txt");
            for(int i = 0; i < userdata.lines.size(); i++) {
                outFile << userdata.lines[i] << endl;
            }
            lines = userdata.lines;
            outFile.close();
        }
        
        for(int i = 0; i < lines.size(); i++) {
            // 选择较长的直线，在该直线上采样
            vector<Point2> line_segments;
            double step = LINES_INTERVAL;
            double dx = lines[i][2] - lines[i][0];
            double dy = lines[i][3] - lines[i][1];
            double length = sqrt(dx * dx + dy * dy);
            double stepX = step * dx / length;
            double stepY = step * dy / length;
            for(double x = lines[i][0], y = lines[i][1]; (x - lines[i][2]) * stepX < 0; x += stepX, y += stepY) {
                line_segments.emplace_back(x, y);
            }
            line_segments.emplace_back(lines[i][2], lines[i][3]);
            selected_lines.emplace_back(line_segments);
        }
    }
    return selected_lines;
}

// 特征点监测
const vector<Point2> & ImageData::getFeaturePoints() const {
    if(feature_points.empty()) {
        FeatureController::detect(getGreyImage(), feature_points, feature_descriptors);
    }
    return feature_points;
}
const vector<FeatureDescriptor> & ImageData::getFeatureDescriptors() const {
    if(feature_descriptors.empty()) {
        FeatureController::detect(getGreyImage(), feature_points, feature_descriptors);
    }
    return feature_descriptors;
}

void ImageData::clear() {
    img.release();
    grey_img.release();
    img_lines.clear();
    feature_points.clear();
    feature_descriptors.clear();
}

// 鼠标事件dz
void onMouse(int event, int x, int y, int flags, void *param) {
    UserData & userdata = *(UserData*)param;
    
    switch (event)
    {
            break;
            //左键按下消息
        case EVENT_LBUTTONDOWN:
        {
            userdata.img.copyTo(userdata.tmp);
            userdata.pre = Point(x, y);
        }
            break;
            //左键抬起消息
        case EVENT_LBUTTONUP:
        {
            userdata.isDelete = true; //置标识符为false
            userdata.cur = Point(x, y);
            line(userdata.img, userdata.pre, userdata.cur, userdata.g_rng.uniform(0, 255), 2, LINE_AA);
            userdata.lines.emplace_back(userdata.pre.x, userdata.pre.y, userdata.cur.x, userdata.cur.y);
        }
            break;
        case EVENT_RBUTTONDOWN:
        {
            if (userdata.isDelete) {
                userdata.tmp.copyTo(userdata.img);
                userdata.lines.pop_back();
                userdata.isDelete = false;
            }
        }
            break;
    }
}