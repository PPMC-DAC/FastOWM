/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "main_options.h"

int setVerbosity(int level)
{
  if((level<vl0) || (level>vl3))
    return 1;

  vFlag = level;
  return 0;
}

void usage(char **argv)
{
	dprintf(vl0,help_text,vl0,vl1,vl2,vl3,vl0,vl1,vl2,vl3,DEFAULT_INPUT_FILE, DEFAULT_NUM_POINTS, DEFAULT_NUM_BLOCK, DEFAULT_DECIMATION, DEFAULT_JUMP_ACEPTABLE, DEFAULT_JUMP_GOOD, DEFAULT_JUMP_PERFECT, DEFAULT_OUTPUT_DIR, DEFAULT_OUTPUT_FILE,DEFAULT_OUTPUT_OBSTACLES_FILE,DEFAULT_NORMALS_RADIUS);
}

void mainOptions_setDefaults(){
    main_options.vValue = DEFAULT_VERBOSITY;
    setVerbosity(main_options.vValue);
    int ret = strlcpy(main_options.inputFileName, DEFAULT_INPUT_FILE,sizeof(main_options.inputFileName));
    if(ret>=sizeof(main_options.inputFileName)) dprintf(vl0,"Input file name was truncated\n");
    main_options.userFileName = 0;
    main_options.userOutputDirName = 0;
    main_options.numPoints  = DEFAULT_NUM_POINTS;
    main_options.num_block = DEFAULT_NUM_BLOCK;
    main_options.decimation = DEFAULT_DECIMATION;
    main_options.jump_aceptable = DEFAULT_JUMP_ACEPTABLE;
    main_options.jump_good = DEFAULT_JUMP_GOOD;
    main_options.jump_perfect = DEFAULT_JUMP_PERFECT;
    main_options.normals_radius = DEFAULT_NORMALS_RADIUS;
    main_options.octree_min_radius = DEFAULT_OCTREE_MIN_RADIUS;
    main_options.do_octree_decimation = DEFAULT_OCTREE_DECIMATE;
    main_options.l_min_points_plane = LANDING2_MINUMUM_POINTS_IN_PLANE;
    main_options.l_stdz_too_high = LANDING2_STDZ_TOO_HIGH;
    main_options.l_steep_slope = LANDING2_STEEP_SLOPE;
    main_options.l_good_slope = LANDING2_GOOD_SLOPE;
    main_options.l_min_skid_circle_r = LANDING2_MINIMUM_SKID_CIRCLE_R;
    main_options.l_good_skid_circle_r = LANDING2_GOOD_SKID_CIRCLE_R;
    main_options.l_min_rotor_circle_r = LANDING2_MINIMUM_ROTOR_CIRCLE_R;
    main_options.l_good_rotor_circle_r = LANDING2_GOOD_ROTOR_CIRCLE_R;
    main_options.l_max_obs_outliers = LANDING2_MAX_OBSTACLES_OUTLIERS;
    main_options.l_height_between_skids_rotor = LANDING2_HEIGTH_BETWEEN_SKIDS_ROTOR;
    main_options.l_height_accept_skids = LANDING2_HEIGTH_ACEPTABLE_SKIDS;
    main_options.l_height_good_skids = LANDING2_HEIGTH_GOOD_SKIDS;
    main_options.l_height_accept_skids_grass_conf = LANDING2_HEIGTH_ACEPTABLE_SKIDS_GRASS_CONF;
    main_options.l_height_good_skids_grass_conf = LANDING2_HEIGTH_GOOD_SKIDS_GRASS_CONF;
    main_options.l_height_accept_skids_plain_conf = LANDING2_HEIGTH_ACEPTABLE_SKIDS_PLAIN_CONF;
    main_options.l_height_good_skids_plain_conf = LANDING2_HEIGTH_GOOD_SKIDS_PLAIN_CONF;
    main_options.l_height_accept_rotor = LANDING2_HEIGTH_ACEPTABLE_ROTOR;
    main_options.l_height_good_rotor = LANDING2_HEIGTH_GOOD_ROTOR;
    main_options.l_pi = LANDING2_PI;
    main_options.l_t_suspicious_angle = LANDING2_SUPICIOUS_ANGLE;
    main_options.l_t_too_rough_angle = LANDING2_TOO_ROUGH_ANGLE;
    main_options.l_t_grass_angle = LANDING2_GRASS_ANGLE;
    main_options.l_t_too_rough_angle_dis = LANDING2_TOO_ROUGH_ANGLE_DIS;
    main_options.l_t_grass_angle_dis = LANDING2_GRASS_ANGLE_DIS;
    main_options.l_t_max_obs_outliers = LANDING2_TERRAIN_MAX_OBSTACLES_OUTLIERS;
    main_options.l_a_aperture_good = LANDING2_APERTURE_GOOD_APROX;
    main_options.l_a_radius = LANDING2_RADIUS_GOOD_APROX;
}

int mainOptions_doGetOpt(int argc, char** argv){
    int c,ret;
    int digit_optind = 0;
    
    /*
     * Copy argv
     */
    // allocate memory and copy strings
    char** new_argv = malloc((argc+1) * sizeof *new_argv);
    int new_argc = argc;
    for(int i = 0; i < argc; ++i)
    {
        size_t length = strlen(argv[i])+1;
        new_argv[i] = malloc(length);
        memcpy(new_argv[i], argv[i], length);
    }
    new_argv[argc] = NULL;
    
    //Keep a copy of optind for getopt_long() to work
    int new_optind = optind;
    
    /*
     * Process Priority arguments, this will be processed first, check for all of them
     */
    while (1) {
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;

        c = getopt_long(new_argc, new_argv, "hDF:f:o:v:n:b:d:j:J:S:N:r:",
                        long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            case 'h' : 
                usage(argv);
                exit(0);
                break;
            case 'F' :
                ret = mainOptions_readFromFile(optarg);
                dprintf(vl0,"Options loaded from file\n");
                break;
            case 'f':
            case 'o':
            case 'v':
            case 'n':
            case 'b':
            case 'd':
            case 'j':
            case 'J':
            case 'S':
            case 'N':
            case 'r':
            case 'D':
            case 1000 :
            case 1001 :
            case 1002 :
            case 1003 :   
            case 1004 :
            case 1005 :
            case 1006 :
            case 1007 :
            case 1008 :
            case 1009 :
            case 1010 :
            case 1011 :
            case 1012 :
            case 1013 :
            case 1014 :
            case 1015 :
            case 1016 :
            case 1017 :
            case 1018 :
            case 1019 :
            case 1020 :
            case 1021 :
            case 1022 :
            case 1023 :
            case 1024 :
            case 1025 :
            case 1026 :
                break;
            case '?':
                break;
            default:
                break;//printf("?? getopt returned character code 0%o ??\n", c);
            }
        }
    
    // free memory for new_argv
    for(int i = 0; i < argc; ++i)
    {
        free(new_argv[i]);
    }
    free(new_argv);
    
    //reset optind
    optind = new_optind;
    
    /*
     * Process other arguments, check for all of them
     */
    
    while (1) {
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;

        c = getopt_long(argc, argv, "hDF:f:o:v:n:b:d:j:J:S:N:r:",
                        long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            case 'h' :
            case 'F' :
                break;
            case 'v':{
       		char *cvValue = optarg;
    	  	main_options.vValue = atoi(cvValue);
                if(setVerbosity(main_options.vValue) != 0){
                    dprintf(vl0,"Verbosity value v (%d) incorrect\n",main_options.vValue);
                    usage(argv);
                    exit(1);
    	  	}
                break;
            }
            case 'o':{
    	  	char *cvValue = optarg;
                ret = strlcpy(main_options.outputDirName, cvValue,sizeof(main_options.outputDirName));
                if(ret>=sizeof(main_options.outputDirName)) dprintf(vl0,"Output dir name was truncated\n");
                dprintf(vl3,"Output directory : %s\n",main_options.outputDirName);
                main_options.userOutputDirName = 1;
                break;
            }
            case 'f':{
       		char *cvValue = optarg;
                ret = strlcpy(main_options.inputFileName, cvValue,sizeof(main_options.inputFileName));
                if(ret>=sizeof(main_options.inputFileName)) dprintf(vl0,"Input file name was truncated\n");
                dprintf(vl3,"Point cloud file to read: %s\n",main_options.inputFileName);
                main_options.userFileName = 1;
                break;
            }
            case 'n':{
       		char *n = optarg;
		main_options.numPoints = (int)atoi(n);
		dprintf(vl3,"Points to read %d\n",main_options.numPoints);
                break;
            }
            case 'b':{
       		char *b = optarg;
                main_options.num_block = (int)atoi(b);
                dprintf(vl3,"Block size set to %d\n",main_options.num_block);
     	    	break;
            }
            case 'd':{
       		char *d = optarg;
                main_options.decimation = (int)atoi(d);
                dprintf(vl3,"Decimation set to %d\n",main_options.decimation);
     	    	break;
            }
            case 'j':{
       		char *j_a = optarg;
 	  	main_options.jump_aceptable = (int)atoi(j_a);
		dprintf(vl3,"Jump on acceptable set to %d\n",main_options.jump_aceptable);
		break;
            }
            case 'J':{
                char *j_g = optarg;
 	  	main_options.jump_good = (int)atoi(j_g);
                dprintf(vl3,"Jump on good set to %d\n",main_options.jump_good);
		break;
            }
            case 'S':{
       		char *j_p = optarg;
 	  	main_options.jump_perfect = (int)atoi(j_p);
                dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
		break;
            }
            case 'N' :
                {
       		char *normR = optarg;
 	  	main_options.normals_radius = (float)atof(normR);
                dprintf(vl3,"Radius to calculate normals set to %.2f\n",main_options.normals_radius);
		break;
            }
            case 'r' :
                {
       		char *minOctR = optarg;
 	  	main_options.octree_min_radius = (float)atof(minOctR);
                dprintf(vl3,"Minimum octree radius set to %.2f\n",main_options.octree_min_radius);
		break;
            }
            case 'D' :
                {
       		//char *minOctR = optarg;
 	  	main_options.do_octree_decimation = 1;//(float)atof(minOctR);
                dprintf(vl3,"YES doing octree decimation\n");
		break;
            }
            case 1000 :
                {
       		char *mpp = optarg;
 	  	main_options.l_min_points_plane = (int)atoi(mpp);
                dprintf(vl3,"Minimum points for a plane set to %d\n",main_options.l_min_points_plane);
		break;
            }
            case 1001 :
                {
       		char *sth = optarg;
 	  	main_options.l_stdz_too_high = (float)atof(sth);
                dprintf(vl3,"Too high stdZ set to %.2f\n",main_options.l_stdz_too_high);
		break;
            }
            case 1002 :
                {
       		char *ss = optarg;
 	  	main_options.l_steep_slope = (float)atof(ss);
                dprintf(vl3,"Steep slope set to %.2f\n",main_options.l_steep_slope);
		break;
            }
            case 1003 :
                {
       		char *gs = optarg;
 	  	main_options.l_good_slope = (float)atof(gs);
                dprintf(vl3,"Good slope set to %.2f\n",main_options.l_good_slope);
		break;
            }
            case 1004 :
                {
       		char *mscr = optarg;
 	  	main_options.l_min_skid_circle_r = (float)atof(mscr);
                dprintf(vl3,"Min skid circle set to %.2f\n",main_options.l_min_skid_circle_r);
		break;
            }
            case 1005 :
                {
       		char *gscr = optarg;
 	  	main_options.l_good_skid_circle_r = (float)atof(gscr);
                dprintf(vl3,"Good skid circle set to %.2f\n",main_options.l_good_skid_circle_r);
		break;
            }
            case 1006 :
                {
       		char *mrcr = optarg;
 	  	main_options.l_min_rotor_circle_r = (float)atof(mrcr);
                dprintf(vl3,"Min rotor circle set to %.2f\n",main_options.l_min_rotor_circle_r);
		break;
            }
            case 1007 :
                {
       		char *grcr = optarg;
 	  	main_options.l_good_rotor_circle_r = (float)atof(grcr);
                dprintf(vl3,"Good rotor circle set to %.2f\n",main_options.l_good_rotor_circle_r);
		break;
            }
            case 1008 :
                {
       		char *hbsr = optarg;
 	  	main_options.l_height_between_skids_rotor = (float)atof(hbsr);
                dprintf(vl3,"Height between skids and rotor set to %.2f\n",main_options.l_height_between_skids_rotor);
		break;
            }
            case 1009 :
                {
       		char *has = optarg;
 	  	main_options.l_height_accept_skids = (float)atof(has);
                dprintf(vl3,"Obstacle height acceptable for skids set to %.2f\n",main_options.l_height_accept_skids);
		break;
            }
            case 1010 :
                {
       		char *hgs = optarg;
 	  	main_options.l_height_good_skids = (float)atof(hgs);
                dprintf(vl3,"Obstacle height good for skids set to %.2f\n",main_options.l_height_good_skids);
		break;
            }
            case 1011 :
                {
       		char *hasgc = optarg;
 	  	main_options.l_height_accept_skids_grass_conf = (float)atof(hasgc);
                dprintf(vl3,"Obstacle height acceptable for skids on grass set to %.2f\n",main_options.l_height_accept_skids_grass_conf);
		break;
            }
            case 1012 :
                {
       		char *hgsgc = optarg;
 	  	main_options.l_height_good_skids_grass_conf = (float)atof(hgsgc);
                dprintf(vl3,"Obstacle height good for skids on grass set to %.2f\n",main_options.l_height_good_skids_grass_conf);
		break;
            }
            case 1013 :
                {
       		char *haspc = optarg;
 	  	main_options.l_height_accept_skids_plain_conf = (float)atof(haspc);
                dprintf(vl3,"Obstacle height acceptable for skids on plain set to %.2f\n",main_options.l_height_accept_skids_plain_conf);
		break;
            }
            case 1014 :
                {
       		char *hgspc = optarg;
 	  	main_options.l_height_good_skids_plain_conf = (float)atof(hgspc);
                dprintf(vl3,"Obstacle height good for skids on plain set to %.2f\n",main_options.l_height_good_skids_plain_conf);
		break;
            }
            case 1015 :
                {
       		char *har = optarg;
 	  	main_options.l_height_accept_rotor = (float)atof(har);
                dprintf(vl3,"Obstacle height acceptable for rotor set to %.2f\n",main_options.l_height_accept_rotor);
		break;
            }
            case 1016 :
                {
       		char *hgr = optarg;
 	  	main_options.l_height_good_rotor = (float)atof(hgr);
                dprintf(vl3,"Obstacle height good for rotor set to %.2f\n",main_options.l_height_good_rotor);
		break;
            }
            case 1017 :
                {
       		char *lpi= optarg;
 	  	main_options.l_pi = (float)atof(lpi);
                dprintf(vl3,"PI set to %.2f\n",main_options.l_pi);
		break;
            }
            case 1018 :
                {
       		char *tsa = optarg;
 	  	main_options.l_t_suspicious_angle = (float)atof(tsa);
                dprintf(vl3,"Suspicious angle for terrain set to: %.2f\n",main_options.l_t_suspicious_angle);
		break;
            }
            case 1019 :
                {
       		char *tsa = optarg;
 	  	main_options.l_t_too_rough_angle = (float)atof(tsa);
                dprintf(vl3,"Too rough angle for terrain (mean) set to: %.2f\n",main_options.l_t_too_rough_angle);
		break;
            }
            case 1020 :
                {
       		char *tga = optarg;
 	  	main_options.l_t_grass_angle = (float)atof(tga);
                dprintf(vl3,"Grass angle for terrain (mean) set to: %.2f\n",main_options.l_t_grass_angle);
		break;
            }
            case 1021 :
                {
       		char *tsad = optarg;
 	  	main_options.l_t_too_rough_angle_dis = (float)atof(tsad);
                dprintf(vl3,"Too rough angle for terrain (dispersion) set to: %.2f\n",main_options.l_t_too_rough_angle_dis);
		break;
            }
            case 1022 :
                {
       		char *tgad = optarg;
 	  	main_options.l_t_grass_angle_dis = (float)atof(tgad);
                dprintf(vl3,"Grass angle for terrain (dispersion) set to: %.2f\n",main_options.l_t_grass_angle_dis);
		break;
            }
            case 1023 :
                {
       		char *aag = optarg;
 	  	main_options.l_a_aperture_good = (float)atof(aag);
                dprintf(vl3,"Good aperture angle for approach set to: %.2f\n",main_options.l_a_aperture_good );
		break;
            }
            case 1024 :
                {
       		char *ar = optarg;
 	  	main_options.l_a_radius = (float)atof(ar);
                dprintf(vl3,"Radius to search for approach obstacles set to: %.2f\n",main_options.l_a_radius );
		break;
            }
            case 1025 :
                {
       		char *tmoo = optarg;
 	  	main_options.l_t_max_obs_outliers = (int)atoi(tmoo);
                dprintf(vl3,"Max obstacles outliers on the terrain set to: %d\n",main_options.l_t_max_obs_outliers );
		break;
            }
            case 1026 :
                {
       		char *moo = optarg;
 	  	main_options.l_max_obs_outliers = (int)atoi(moo);
                dprintf(vl3,"Max obstacles outliers on the landing zone set to: %d\n",main_options.l_max_obs_outliers );
		break;
            }
            case '?':
                {
                if (optopt)  printf("Unrecognized short option '%c'\n", optopt);
                else  printf("Unrecognized long option \"%s\"\n", argv[this_option_optind]);
                break;
            }
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
            }
    }
    
    return optind;   
}


void mainOptions_print(int verbosity_level){
    dprintf(verbosity_level,"Verbosity set to: %d\n",main_options.vValue);
    dprintf(verbosity_level,"Input Filename set to: %s\n",main_options.inputFileName);
    dprintf(verbosity_level,"Output Directory name set to: %s\n",main_options.outputDirName);
    dprintf(verbosity_level,"Points to read: %d\n",main_options.numPoints);
    dprintf(verbosity_level,"Block size set to: %d\n",main_options.num_block);
    dprintf(verbosity_level,"Decimation set to: %d\n",main_options.decimation);
    dprintf(verbosity_level,"Jump on acceptable set to: %d\n",main_options.jump_aceptable);
    dprintf(verbosity_level,"Jump on good set to: %d\n",main_options.jump_good);
    dprintf(verbosity_level,"Jump on perfect set to: %d\n",main_options.jump_perfect); 
    dprintf(verbosity_level,"Radius to calculate normals set to: %.2f\n",main_options.normals_radius);
    dprintf(verbosity_level,"Octree minimum radius set to: %.2f\n",main_options.octree_min_radius);
    if(main_options.do_octree_decimation){
        dprintf(verbosity_level,"YES doing octree decimation\n");
    }else{
        dprintf(verbosity_level,"NO doing octree decimation\n");
    }
    dprintf(verbosity_level,"Minimum number of points to build a plane: %d\n",main_options.l_min_points_plane);
    dprintf(verbosity_level,"Standard deviation of a plane in Z considered too high: %.2f\n", main_options.l_stdz_too_high);
    dprintf(verbosity_level,"Steep slope limit: %.2f\n",main_options.l_steep_slope);
    dprintf(verbosity_level,"Good slope limit: %.2f\n",main_options.l_good_slope);
    dprintf(verbosity_level,"Minimum skid radius: %.2f\n",main_options.l_min_skid_circle_r);
    dprintf(verbosity_level,"Good skid radius: %.2f\n",main_options.l_good_skid_circle_r);
    dprintf(verbosity_level,"Minimum rotor radius: %.2f\n",main_options.l_min_rotor_circle_r);
    dprintf(verbosity_level,"Good rotor radius: %.2f\n",main_options.l_good_rotor_circle_r);
    dprintf(verbosity_level,"Max obstacles outliers on the landing zone: %d\n",main_options.l_max_obs_outliers);
    dprintf(verbosity_level,"Height between rotor and skids: %.2f\n",main_options.l_height_between_skids_rotor);
    dprintf(verbosity_level,"Acceptable obstacle heigh for the skids: %.2f\n",main_options.l_height_accept_skids);
    dprintf(verbosity_level,"Good obstacle heigh for the skids: %.2f\n",main_options.l_height_good_skids);
    dprintf(verbosity_level,"Acceptable obstacle heigh for the skids on grass: %.2f\n",main_options.l_height_accept_skids_grass_conf);
    dprintf(verbosity_level,"Good obstacle heigh for the skids on grass: %.2f\n",main_options.l_height_good_skids_grass_conf);
    dprintf(verbosity_level,"Acceptable obstacle heigh for the skids on plain: %.2f\n",main_options.l_height_accept_skids_plain_conf);
    dprintf(verbosity_level,"Good obstacle heigh for the skids on plain: %.2f\n",main_options.l_height_good_skids_plain_conf);
    dprintf(verbosity_level,"Acceptable obstacle heigh for the rotor: %.2f\n",main_options.l_height_accept_rotor);
    dprintf(verbosity_level,"Good obstacle heigh for the rotor: %.2f\n",main_options.l_height_good_rotor);
    dprintf(verbosity_level,"Value of PI: %.2f\n",main_options.l_pi);
    dprintf(verbosity_level,"Suspicious angle for obstacles: %.2f\n",main_options.l_t_suspicious_angle);
    dprintf(verbosity_level,"Angle considered too rough for terrain (mean): %.2f\n",main_options.l_t_too_rough_angle);
    dprintf(verbosity_level,"Angle considered grass for terrain (mean): %.2f\n",main_options.l_t_grass_angle);
    dprintf(verbosity_level,"Angle considered too rough for terrain (dispersion): %.2f\n",main_options.l_t_too_rough_angle_dis);
    dprintf(verbosity_level,"Angle considered grass for terrain (dispersion): %.2f\n",main_options.l_t_grass_angle_dis);
    dprintf(verbosity_level,"Max obstacles outliers on the terrain: %d\n",main_options.l_t_max_obs_outliers);
    dprintf(verbosity_level,"Aperture angle for a clear approach: %.2f\n",main_options.l_a_aperture_good);
    dprintf(verbosity_level,"Radius to search for approach obstacles: %.2f\n",main_options.l_a_radius);
}


int mainOptions_readFromFile(char * file){
    char variable[MAX_FILE_NAME_LENGTH],value[MAX_FILE_NAME_LENGTH];
    int ret;
    
    FILE * f = fopen(file, "r");
    if(f==NULL){
        dprintf(vl3,"Options file could not be opened : %s\n",file);
        return 0;
    }
    while(!feof(f)) {
        if(fscanf(f, "%s %s%*[^\n]", variable, value) == 2){
           if(strcmp(variable, "verbosity") == 0){
                main_options.vValue = atoi(value);
                if(setVerbosity(main_options.vValue) != 0){
                    printf("Verbosity value v (%d) incorrect\n",main_options.vValue);
                    exit(1);
                }
           } else if(strcmp(variable, "output-dir") == 0) {
               ret = strlcpy(main_options.outputDirName, value,sizeof(main_options.outputDirName));
               if(ret>=sizeof(main_options.outputDirName)) dprintf(vl0,"Output dir name was truncated\n");
               //dprintf(vl3,"Output directory : %s\n",main_options.outputDirName);
               main_options.userOutputDirName = 1;
           } else if(strcmp(variable, "point-cloud-file") == 0) {
                ret = strlcpy(main_options.inputFileName, value,sizeof(main_options.inputFileName));
                if(ret>=sizeof(main_options.inputFileName)) dprintf(vl0,"Input file name was truncated\n");
                //dprintf(vl3,"Point cloud file to read: %s\n",main_options.inputFileName);
                main_options.userFileName = 1;
            }else if(strcmp(variable, "num-points") == 0) {
                main_options.numPoints = (int)atoi(value);
                //dprintf(vl3,"Points to read %d\n",main_options.numPoints);
            }else if(strcmp(variable, "block-size") == 0) {
                main_options.num_block = (int)atoi(value);
                //dprintf(vl3,"Block size set to %d\n",main_options.num_block);
            }else if(strcmp(variable, "decimation") == 0) {
                main_options.decimation = (int)atoi(value);
                //dprintf(vl3,"Decimation set to %d\n",main_options.decimation);
            }else if(strcmp(variable, "jump-acceptable") == 0) {
                main_options.jump_aceptable = (int)atoi(value);
                //dprintf(vl3,"Jump on acceptable set to %d\n",main_options.jump_aceptable);
            }else if(strcmp(variable, "jump-good") == 0) {
                main_options.jump_good = (int)atoi(value);
                //dprintf(vl3,"Jump on good set to %d\n",main_options.jump_good);
            }else if(strcmp(variable, "jump-perfect") == 0) {
                main_options.jump_perfect = (int)atoi(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "normals-radius") == 0) {
                main_options.normals_radius = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "octree-min-radius") == 0) {
                main_options.octree_min_radius = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "do-octree-decimation") == 0) {
                main_options.do_octree_decimation = (int)atoi(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-min-points-plane") == 0) {
                main_options.l_min_points_plane = (int)atoi(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-stdz-too-high") == 0) {
                main_options.l_stdz_too_high = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-good-slope") == 0) {
                main_options.l_good_slope = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-steep-slope") == 0) {
                main_options.l_steep_slope = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-min-skid-circle-r") == 0) {
                main_options.l_min_skid_circle_r = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-good-skid-circle-r") == 0) {
                main_options.l_good_skid_circle_r = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-min-rotor-circle-r") == 0) {
                main_options.l_min_rotor_circle_r= (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-good-rotor-circle-r") == 0) {
                main_options.l_good_rotor_circle_r = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-max-obs-outliers") == 0) {
                main_options.l_max_obs_outliers = (int)atoi(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-between-skids-rotor") == 0) {
                main_options.l_height_between_skids_rotor = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-accept-skids") == 0) {
                main_options.l_height_accept_skids = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-good-skids") == 0) {
                main_options.l_height_good_skids = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-accept-skids-grass-conf") == 0) {
                main_options.l_height_accept_skids_grass_conf = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-good-skids-grass-conf") == 0) {
                main_options.l_height_good_skids_grass_conf = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-accept-skids-plain-conf") == 0) {
                main_options.l_height_accept_skids_plain_conf = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-good-skids-plain-conf") == 0) {
                main_options.l_height_good_skids_plain_conf = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-accept-rotor") == 0) {
                main_options.l_height_accept_rotor = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-height-good-rotor") == 0) {
                main_options.l_height_good_rotor = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-pi") == 0) {
                main_options.l_pi = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-suspicious-angle") == 0) {
                main_options.l_t_suspicious_angle= (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-too-rough-angle") == 0) {
                main_options.l_t_too_rough_angle = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-grass-angle") == 0) {
                main_options.l_t_grass_angle = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-too-rough-angle-dis") == 0) {
                main_options.l_t_too_rough_angle_dis = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-grass-angle-dis") == 0) {
                main_options.l_t_grass_angle_dis = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-t-max-obs-outliers") == 0) {
                main_options.l_t_max_obs_outliers = (int)atoi(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-a-aperture-good") == 0) {
                main_options.l_a_aperture_good = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }else if(strcmp(variable, "l-a-radius") == 0) {
                main_options.l_a_radius = (float)atof(value);
                //dprintf(vl3,"Jump on perfect set to %d\n",main_options.jump_perfect);
            }
        }
        else{
            printf("Error in %s %s\n", variable, value);
        }
    }
    
    fclose(f);
    return 1;
}

