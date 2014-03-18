//
// Copyright 2013-2014, Bradley A. Grantham
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//      http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <ctime>
#include <cerrno>
#include <sys/time.h>
#include <vector>
#include <limits>

#include "world.h"

world *gWorld;

float angle_around = 0.0f;
float angle_over = 0.0f;
float zoom = 0.0f;

unsigned int scene_data_row_pad = 512;
scene_shader_data scene_data;

void load_scene_data(world *w)
{
    unsigned int spheresqrt = (unsigned int)ceil(sqrt(w->sphere_count));
    scene_data_row_pad = 1;
    while(scene_data_row_pad < spheresqrt)
        scene_data_row_pad *= 2;

    get_shader_data(w, scene_data, scene_data_row_pad);
}

void usage(char *progname)
{
    fprintf(stderr, "usage: %s <inputfilename>\n", progname);
    fprintf(stderr, "A PPM file is written as output.");
}

void update_world_camera(world *world, float around, float over, float zoom)
{
    vec3 offset;

    offset.x = zoom * cos(over) * sin(around);
    offset.y = zoom * sin(over);
    offset.z = -zoom * cos(over) * cos(around);

    world->cam.eye = vec3_add(world->scene_center, offset);
    world->cam.yaw = -around;
    world->cam.pitch = - over;
}

int digits(float f)
{
    if(f == trunc(f))
        return 0;
    if(f * 10 == trunc(f * 10))
        return 1;
    if(f * 100 == trunc(f * 100))
        return 2;
    if(f * 1000 == trunc(f * 1000))
        return 3;
    if(f * 10000 == trunc(f * 10000))
        return 4;
    if(f * 100000 == trunc(f * 100000))
        return 5;
    return 6;
}

void dump_array(const std::string field_name, int columns, int rows, float *data, int size, bool rewrite_terminators = false)
{
    printf("scene_data.%s = new Float32Array([\n", field_name.c_str());
    for(int i = 0; i < rows; i++)
        for(int j = 0; j < columns; j++) {
            int which = i * columns + j;
            for(int k = 0; k < size; k++) {
                float f = data[which * size + k];
                if(rewrite_terminators && (f == float(0x7fffffff)))
                    f = float(0xffff);
                if(f == 0.0)
                    printf("0,");
                else
                    printf("%.*f,", digits(f), f);
            }
            if(size == 2)
                printf("0,");
            //printf("\n");
        }
        printf("])\n");
}

int main(int argc, char *argv[])
{
    if(argc < 2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if((!strcmp(argv[1], "-h")) || (!strcmp(argv[1], "--help"))) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if((gWorld = load_world(argv[1])) == NULL) {
        fprintf(stderr, "Cannot set up world.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "loaded\n");
    gWorld->xsub = gWorld->ysub = 1;
    gWorld->cam.fov = to_radians(40.0);
    zoom = gWorld->scene_extent / 2 / sinf(gWorld->cam.fov / 2);
    update_world_camera(gWorld, angle_around, angle_over, zoom);

    load_scene_data(gWorld);

    printf("scene_data_row_pad = %d\n", scene_data_row_pad);

    printf("scene_data = {}\n");
    printf("scene_data.row_pad = %d\n", scene_data_row_pad);
    printf("scene_data.sphere_count = %d\n", scene_data.sphere_count);
    printf("scene_data.sphere_data_rows = %d\n", scene_data.sphere_data_rows);
    printf("scene_data.group_data_rows = %d\n", scene_data.group_data_rows);
    printf("scene_data.group_count = %d\n", scene_data.group_count);
    printf("scene_data.tree_root = %d\n", scene_data.tree_root);

    dump_array("sphere_geometries", scene_data_row_pad, scene_data.sphere_data_rows, scene_data.sphere_geometries, 4);
    dump_array("sphere_colors", scene_data_row_pad, scene_data.sphere_data_rows, scene_data.sphere_colors, 3);

    dump_array("group_boxmin", scene_data_row_pad, scene_data.group_data_rows, scene_data.group_boxmin, 3);
    dump_array("group_boxmax", scene_data_row_pad, scene_data.group_data_rows, scene_data.group_boxmax, 3);
    dump_array("group_children", scene_data_row_pad, scene_data.group_data_rows, scene_data.group_children, 2, true);
    dump_array("group_objects", scene_data_row_pad, scene_data.group_data_rows, scene_data.group_objects, 2, true);
    dump_array("group_hitmiss", scene_data_row_pad, scene_data.group_data_rows, scene_data.group_hitmiss, 2, true);

    printf("scene_data.scene_center = [%f, %f, %f]\n", gWorld->scene_center.x, gWorld->scene_center.y, gWorld->scene_center.z);
    printf("scene_data.scene_extent = %f\n", gWorld->scene_extent);
    printf("scene_data.eye = [%f, %f, %f]\n", gWorld->cam.eye.x, gWorld->cam.eye.y, gWorld->cam.eye.z);
    printf("scene_data.yaw = %f\n", gWorld->cam.yaw);
    printf("scene_data.pitch = %f\n", gWorld->cam.pitch);
    printf("scene_data.roll = %f\n", gWorld->cam.roll);
    printf("scene_data.fov = %f\n", gWorld->cam.fov);
    printf("scene_data.background = [%f, %f, %f]\n", gWorld->background.x, gWorld->background.y, gWorld->background.z);
    printf("window.scene_data = scene_data\n");

    exit(EXIT_SUCCESS);
}
