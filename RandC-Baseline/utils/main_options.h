/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   option_file.h
 * Author: oscar.garcia
 *
 * Created on November 15, 2018, 11:35 AM
 */

#ifndef MAIN_OPTIONS_H

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <bsd/string.h> //in Ubuntu this is needed for strlcpy

#define MAIN_OPTIONS_H

#ifdef __cplusplus
extern "C" {
#endif



/*
 * Debugging (verbosity level)
 */
#define vl0 0
#define vl1 1
#define vl2 2
#define vl3 3
int vFlag;
#define DEFAULT_VERBOSITY vl0


int setVerbosity(int level);

/** dprintf: all messages whose priority level is equal
 * or higher than vFlag (set by setVerbosity()) will
 * be printed
 */

#define __USE_DPRINTF  //Comment to remove any output on screen

#ifdef __USE_DPRINTF
#define dprintf(level,...) \
{\
   if(vFlag >= level) printf(__VA_ARGS__);\
}
#else
#define dprintf(level,...) " "
#endif


/**
 * Define defaults
 */
#define MAX_FILE_NAME_LENGTH 254

#define DEFAULT_NORMALS_RADIUS 1.5
#define DEFAULT_OCTREE_MIN_RADIUS 1.5
#define DEFAULT_OCTREE_DECIMATE 0

#define DEFAULT_INPUT_FILE "input.xyz"
#define DEFAULT_OUTPUT_DIR ".result.landing/" //Designed to include point (.) and slash (/) to be concatenated between input file path and output files names
#define DEFAULT_OUTPUT_FILE "results2.xyz"
#define DEFAULT_OUTPUT_OBSTACLES_FILE "obstacles.txt"
#define DEFAULT_NUM_POINTS 100000
#define DEFAULT_NUM_BLOCK 64
#define DEFAULT_DECIMATION 1
#define DEFAULT_JUMP_ACEPTABLE 8
#define DEFAULT_JUMP_GOOD 32
#define DEFAULT_JUMP_PERFECT 256

#define LANDING2_MINUMUM_POINTS_IN_PLANE 5
#define LANDING2_STDZ_TOO_HIGH 1

#define LANDING2_STEEP_SLOPE 0.1 //angle in degrees = atan(angle in percent / 100%)
#define LANDING2_GOOD_SLOPE 0.08 //angle in degrees = atan(angle in percent / 100%)

#define LANDING2_MINIMUM_SKID_CIRCLE_R 2.0
#define LANDING2_GOOD_SKID_CIRCLE_R 11.0
#define LANDING2_MINIMUM_ROTOR_CIRCLE_R 11.0
#define LANDING2_GOOD_ROTOR_CIRCLE_R 12.5
#define LANDING2_HEIGTH_BETWEEN_SKIDS_ROTOR 4.0
#define LANDING2_MAX_OBSTACLES_OUTLIERS 4.0

#define LANDING2_HEIGTH_ACEPTABLE_SKIDS 0.5
#define LANDING2_HEIGTH_GOOD_SKIDS 0.3
#define LANDING2_HEIGTH_ACEPTABLE_SKIDS_GRASS_CONF 0.5
#define LANDING2_HEIGTH_GOOD_SKIDS_GRASS_CONF 0.4
#define LANDING2_HEIGTH_ACEPTABLE_SKIDS_PLAIN_CONF 0.4
#define LANDING2_HEIGTH_GOOD_SKIDS_PLAIN_CONF 0.2
#define LANDING2_HEIGTH_ACEPTABLE_ROTOR 1.5
#define LANDING2_HEIGTH_GOOD_ROTOR 0.5


#define LANDING2_PI 3.1416
#define LANDING2_SUPICIOUS_ANGLE 0.7854 //quarter pi, 45º
#define LANDING2_TOO_ROUGH_ANGLE 3.1416 //pi, 180º
#define LANDING2_GRASS_ANGLE 1.5708 //half pi, 90º
//#define LANDING2_GRASS_ANGLE 1,7453 // 100º
#define LANDING2_TOO_ROUGH_ANGLE_DIS 3.1416 //pi, 180º
#define LANDING2_GRASS_ANGLE_DIS 1.5708 //half pi, 90º
#define LANDING2_TERRAIN_MAX_OBSTACLES_OUTLIERS 0

#define LANDING2_APERTURE_GOOD_APROX 160.0
#define LANDING2_RADIUS_GOOD_APROX 30.0

/*
 * Define options (best for reading from file)
 */

    struct main_options{
        int vValue;
        char inputFileName[MAX_FILE_NAME_LENGTH];
        char userFileName;
        char outputDirName[MAX_FILE_NAME_LENGTH];
        char userOutputDirName;
        int numPoints;
        int num_block;
        int decimation;
        int jump_aceptable;
        int jump_good;
        int jump_perfect;
        float normals_radius;
        float octree_min_radius;
        int do_octree_decimation;
        int l_min_points_plane;
        float l_stdz_too_high;
        float l_steep_slope;
        float l_good_slope;
        float l_min_skid_circle_r;
        float l_good_skid_circle_r;
        float l_min_rotor_circle_r;
        float l_good_rotor_circle_r;
        int l_max_obs_outliers;
        float l_height_between_skids_rotor;
        float l_height_accept_skids;
        float l_height_good_skids;
        float l_height_accept_skids_grass_conf;
        float l_height_good_skids_grass_conf;
        float l_height_accept_skids_plain_conf;
        float l_height_good_skids_plain_conf;
        float l_height_accept_rotor;
        float l_height_good_rotor;
        float l_pi;
        float l_t_suspicious_angle;
        float l_t_too_rough_angle;
        float l_t_grass_angle;
        float l_t_too_rough_angle_dis;
        float l_t_grass_angle_dis;
        int l_t_max_obs_outliers;
        float l_a_aperture_good;
        float l_a_radius;
    }main_options;

    static struct option long_options[] ={
                    {"help", no_argument, 0, 'h'},
                    {"verbosity", required_argument, 0, 'v'},
                    {"options-file", required_argument, 0, 'F'},
                    {"output-dir",  required_argument, 0, 'o'},
                    {"point-cloud-file",  required_argument, 0, 'f'},
                    {"num-points",  required_argument, 0, 'n'},
                    {"block-size",    required_argument, 0, 'b'},
                    {"decimation",    required_argument, 0, 'd'},
                    {"jump-acceptable",    required_argument, 0, 'j'},
                    {"jump-good",    required_argument, 0, 'J'},
                    {"jump-perfect",    required_argument, 0, 'S'},
                    {"normals-radius", required_argument, 0, 'N'},
                    {"octree-min-radius", required_argument, 0, 'r'},
                    {"do-octree-decimation", required_argument, 0, 'D'},
                    {"l-min-points-plane", required_argument, 0, 1000},
                    {"l-stdz-too-high", required_argument, 0, 1001},
                    {"l-steep-slope", required_argument, 0, 1002},
                    {"l-good-slope", required_argument, 0, 1003},
                    {"l-min-skid-circle-r", required_argument, 0, 1004},
                    {"l-good-skid-circle-r", required_argument, 0, 1005},
                    {"l-min-rotor-circle-r", required_argument, 0, 1006},
                    {"l-good-rotor-circle-r", required_argument, 0, 1007},
                    {"l-max-obs-outliers", required_argument, 0, 1026},
                    {"l-height-between-skids-rotor", required_argument, 0, 1008},
                    {"l-height-accept-skids", required_argument, 0, 1009},
                    {"l-height-good-skids", required_argument, 0, 1010},
                    {"l-height-accept-skids-grass-conf", required_argument, 0, 1011},
                    {"l-height-good-skids-grass-conf", required_argument, 0, 1012},
                    {"l-height-accept-skids-plain-conf", required_argument, 0, 1013},
                    {"l-height-good-skids-plain-conf", required_argument, 0, 1014},
                    {"l-height-accept-rotor", required_argument, 0, 1015},
                    {"l-height-good-rotor", required_argument, 0, 1016},
                    {"l-pi", required_argument, 0, 1017},
                    {"l-t-suspicious-angle", required_argument, 0, 1018},
                    {"l-t-too-rough-angle", required_argument, 0, 1019},
                    {"l-t-grass-angle", required_argument, 0, 1020},
                    {"l-t-too-rough-angle-dis", required_argument, 0, 1021},
                    {"l-t-grass-angle-dis", required_argument, 0, 1022},
                    {"l-t-max-obs-outliers", required_argument, 0, 1025},
                    {"l-a-aperture-good", required_argument, 0, 1023},
                    {"l-a-radius", required_argument, 0, 1024},
                    {0,         0,                 0,  0 }
               };

    static char * help_text = "This program classifies LiDAR points according how adequate they are as a center for a helicopter landing site.\n"
                        "As an input it requires a LiDAR point cloud file, in text format with at least 4 fields, for X Y Z in double point and Intensity in integer.\n"
                        "It was tested in an Ubuntu Linux.\n"
                        "Its options are:\n"
                        "\t-h, --help \t\t\tThis message\n"
                        "\t-F, --options-file=[+STRING] \t\t\tAn options file where the parameters of the program are read from, command line paraters will have priority if used\n"
                        "\t-v, --verbosity=[+NUM] \t\t\tThe verbose level, %d less information, %d recommended, %d debug, %d all\n"
                        "\t-o, --output_dir=[+STRING] \t\t\tThe output directory (a / will be added to the name at the end), if this option is not used the default behaviour is commented later in this document\n"
                        "\t-f, --point_cloud_file=[+STRING] \t\t\tThe input file, the default is %s\n"
                        "\t-n, --num_points=[+NUM] \t\t\tThe number of points to be read and classified, if the file ends before all these points are read it is OK, but if it is not especified the default is %d\n"
                        "\t-b, --block_size=[+NUM] \t\t\tThe block sixe for the omp worh share, the default is %d\n"
                        "\t-d, --decimation=[+NUM] \t\t\tThe decimation of the input file, for example a decimation of 5 will pick 1 in 5 points to be classified, the default is %d, no decimation\n"
                        "\t-j, --jump_acceptable=[+NUM] \t\t\tAfter an acceptable landing site is found, the number of points in its neighbourhood that will not be searched for a new landing site, the default is %d\n"
                        "\t-J, --jump_good=[+NUM] \t\t\tAfter a good landing site is found, the number of points in its neighbourhood that will not be searched for a new landing site, the default is %d\n"
                        "\t-S, --jump_perfect=[+NUM] \t\t\tAfter a perfect landing site is found, the number of points in its neighbourhood that will not be searched for a new landing site, the default is %d\n"
                        "\t-N, --normals-radius=[+NUM] \t\t\tThe radius that will be used to calculate the normals, the default is %.2f. For each point its normal is calculated using its neighbours; the larder the point cloud density the larger this radius should be, to increase precision and speed\n"
                        "\t-r, --octree-min-radius=[+NUM] \t\t\tThe minimum radius of the octree, the point where it no longer divides, the maximum radius of its leafs, the default is %.2f.  Because the octree halves its radius at each level the actual radius may not be the same. The larger the point density of the cloud, the smaller this radius may be\n"
                        "\t-D, --do-octree-decimation \t\t\tDecimate according to the octree, keeping the highest and lowest point of each leaf. The larger the minumum radius of the octree more points will be decimated. This decimation means that all the points will be tested as landing zones, but the points where obstacles are looked for are reduced. Because the reduction is done in the octree, its is spatially balanced and should not provoke many errors as long as the octree radius is well selected\n"
                        "\t--l-min-points-plane=[+NUM] \t\t\tHow many points are needed in a neighbourhood to fit a plane\n"
                        "\t--l-stdz-too-high=[+NUM] \t\t\tNOT CURRENTLY USED The maximum standard deviation on Z axis of the points fitting a plane, higher would mean the plane is too rough for landing\n"
                        "\t--l-steep-slope=[+NUM] \t\t\tThe slope from which a landing is possible but difficult, in angle in radians (angle in degrees = atan(angle in percent / 100\%))\n"
                        "\t--l-good-slope=[+NUM] \t\t\tThe slope from which a landing is no longer possible, in angle in radians (angle in degrees = atan(angle in percent / 100\%))\n"
                        "\t--l-min-skid-circle-r=[+NUM] \t\t\tThe radius of the minimum circle for landing the skids (in meters)\n"
                        "\t--l-good-skid-circle-r=[+NUM] \t\t\tThe radius of the good circle for landing the skids (in meters)\n"
                        "\t--l-min-rotor-circle-r=[+NUM] \t\t\tNOT CURRENTLY USED The radius of the minimum circle for landing the rotor and fuselage (in meters)\n"
                        "\t--l-good-rotor-circle-r=[+NUM] \t\t\tThe radius of the good circle for landing the rotor and fuselage (in meters)\n"
                        "\t--l-max-obs-outliers=[+NUM] \t\t\tThe number of points detected as obstacles in the landing area that are ignored as real obstacles\n"
                        "\t--l-height-between-skids-rotor=[+NUM] \t\t\tNOT CURRENTLY USED The height between the skids an the rotors (in meters)\n"
                        "\t--l-height-accept-skids=[+NUM] \t\t\tThe maximum height of vegetation or obstacles for an acceptable landing of the skids (in meters)\n"
                        "\t--l-height-good-skids=[+NUM] \t\t\tThe maximum height of vegetation or obstacles for a good landing of the skids (in meters)\n"
                        "\t--l-height-accept-skids-grass-conf=[+NUM] \t\t\tThe maximum height of vegetation or obstacles in a grassy zone for an acceptable landing of the skids if the vegetation  or obstacles is though to be higher that the data shows (in meters)\n"
                        "\t--l-height-good-skids-grass-conf=[+NUM] \t\t\tThe maximum height of vegetation or obstacles in a grassy zone for an good landing of the skids if the vegetation  or obstacles is though to be higher that the data shows (in meters)\n"
                        "\t--l-height-accept-skids-plain-conf=[+NUM] \t\t\tThe maximum height of vegetation or obstacles in a plain zone for an acceptable landing of the skids if the vegetation  or obstacles is though to be higher that the data shows (in meters)\n"
                        "\t--l-height-good-skids-plain-conf=[+NUM] \t\t\tThe maximum height of vegetation or obstacles in a plain zone for a good landing of the skids if the vegetation  or obstacles is though to be higher that the data shows (in meters)\n"
                        "\t--l-height-accept-rotor=[+NUM] \t\t\tThe maximum height of vegetation or obstacles for an acceptable landing of the rotors (in meters\n"
                        "\t--l-height-good-rotor=[+NUM] \t\t\tThe maximum height of vegetation or obstacles for a good landing of the rotors (in meters)\n"
                        "\t--l-pi=[+NUM] \t\t\tthe value of PI to be used on computations\n"
                        "\t--l-t-suspicious-angle=[+NUM] \t\t\tThe angle of the normal of an obstacle point where that obstacle is though of being actually higher that the data shows (in radians)\n"
                        "\t--l-t-too-rough-angle=[+NUM] \t\t\tNOT CURRENTLY USED The angle of the normal of each point in a terrain plane to be considered too rough (accumulated as their mean) (in radians)\n"
                        "\t--l-t-grass-angle=[+NUM] \t\t\tNOT CURRENTLY USED The angle of the normal of each point in a terrain plane to be considered grass (accumulated as their mean) (in radians)\n"
                        "\t--l-t-too-rough-angle-dis=[+NUM] \t\t\tThe angle of the normal of each point in a terrain plane to be considered too rough (accumulated as dispersion) (in radians)\n"
                        "\t--l-t-grass-angle-dis=[+NUM] \t\t\tThe angle of the normal of each point in a terrain plane to be considered grass (accumulated as dispersion) (in radians)\n"
                        "\t--l-t-max-obs-outliers=[+NUM] \t\t\tThe number of points detected as obstacles in the terrain that are ignored as real obstacles\n"
                        "\t--l-a-aperture-good=[+NUM] \t\t\tThe aperture of the cone to be used to calculate approach obstacles (the approach angle should be subtracted twice to 180º to get the aperture) (in degrees)\n"
                        "\t--l-a-radius=[+NUM] \t\t\tThe radius of the cone to be used to calculate approach obstacles (in meters)\n"
                        "To change the number of OpenMP threads change the variable OMP_NUM_THREADS in bash.\n"
                        "The output will be saved in a folder with a composite name beginning with the input file name and ending with %s.\n"
                        "The output files are:\n"
                        "\t%s : A point cloud file with X Y Z I Classification\n"
                        "\t%s : A file with the obstacles for each approximation, with a point followed of several pairs of angles in format X Y Z [angle_a angle_b]* ... [angle_a angle_b]*\n"
                        "The values that define the landing zones and their classification must be changed in the source code by changing several #define\n"
                    ;


    void usage(char **argv);
    void mainOptions_setDefaults();
    int mainOptions_doGetOpt(int argc, char** argv);
    void mainOptions_print(int verbosity_level);
    int mainOptions_readFromFile(char * file);




#ifdef __cplusplus
}
#endif

#endif /* MAIN_OPTIONS_H */
