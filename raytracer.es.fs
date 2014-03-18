#version 100

precision highp float;

uniform int which;
uniform int tree_root;
uniform highp vec4 light_dir;
uniform highp vec3 eye;
uniform highp mat4 camera_matrix;
uniform highp mat4 camera_normal_matrix;
uniform highp mat4 object_matrix;
uniform highp mat4 object_inverse;
uniform highp mat4 object_normal_matrix;
uniform highp mat4 object_normal_inverse;
uniform float aspect;
uniform float fov;
uniform highp vec3 background_color;
uniform sampler2D sphere_geometry;
uniform int sphere_data_row_size;
uniform int sphere_data_rows;
uniform sampler2D sphere_colors;
uniform int group_data_row_size;
uniform int group_data_rows;
uniform sampler2D group_boxmin;
uniform sampler2D group_boxmax;
uniform sampler2D group_children;
uniform sampler2D group_objects;
uniform sampler2D group_hitmiss;

varying highp vec2 ftex;

struct ray {
    highp vec4 o;
    highp vec4 d;
};

struct surface_hit
{
    float t;
    int which;
    vec3 color;
};

const float infinitely_far = 10000000.0;

surface_hit surface_hit_init(vec3 background_color)
{
    return surface_hit(infinitely_far, -1, background_color);
}

void set_bad_hit(inout surface_hit hit, float r, float g, float b)
{
    hit.color = vec3(r, g, b);
    hit.t = 1.0;
}

struct range
{
    float t0, t1;
};

range range_full()
{
    return range(-100000000.0, 100000000.0);
}

range range_intersect(in range r1, in range r2)
{
    float t0 = max(r1.t0, r2.t0);
    float t1 = min(r1.t1, r2.t1);
    return range(t0, t1);
}

bool range_is_empty(in range r)
{
    return r.t0 >= r.t1;
}

const float NO_t = -1.0;

range box_intersect(in highp vec3 boxmin, in highp vec3 boxmax, in ray theray)
{
    // XXX w ignored

    range r = range_full();

    float t0, t1;

    t0 = (boxmin.x - theray.o.x) / theray.d.x;
    t1 = (boxmax.x - theray.o.x) / theray.d.x;
    if(theray.d.x >= 0.0)
        r = range_intersect(r, range(t0, t1));
    else
        r = range_intersect(r, range(t1, t0));

    t0 = (boxmin.y - theray.o.y) / theray.d.y;
    t1 = (boxmax.y - theray.o.y) / theray.d.y;
    if(theray.d.y >= 0.0)
        r = range_intersect(r, range(t0, t1));
    else
        r = range_intersect(r, range(t1, t0));

    t0 = (boxmin.z - theray.o.z) / theray.d.z;
    t1 = (boxmax.z - theray.o.z) / theray.d.z;
    if(theray.d.z >= 0.0)
        r = range_intersect(r, range(t0, t1));
    else
        r = range_intersect(r, range(t1, t0));

    return r;
}

range sphere_intersect(in highp vec3 center, in float radius, in ray theray)
{
    // XXX w ignored
    highp vec3 diff = theray.o.xyz - center;

    float b = 2.0 * dot(theray.d.xyz, diff);

    float radicand = b * b - 4.0 * (dot(diff, diff) - radius * radius);

    if(radicand < 0.0)
        return range(100000000.0, -100000000.0);

    float t0 = (-b - sqrt(radicand)) / 2.0;
    float t1 = (-b + sqrt(radicand)) / 2.0;
    return range(t0, t1);
}

float sphere_intersect(in highp vec3 center, in float radius, in ray theray, in range r)
{
    // XXX w ignored
    highp vec3 diff = theray.o.xyz - center;

    float b = 2.0 * dot(theray.d.xyz, diff);

    float radicand = b * b - 4.0 * (dot(diff, diff) - radius * radius);
    if(radicand < 0.0)
        return NO_t;

    float t1 = (-b + sqrt(radicand)) / 2.0;

    if(t1 < r.t0)
         return NO_t;

    float t0 = (-b - sqrt(radicand)) / 2.0;

    if(t0 > r.t1)
        return NO_t;

    if(t0 < r.t0)
        return t1;

    return t0;
}

struct group {
    bool is_branch;
    int g1, g2;
    int start, count;
    highp vec3 boxmin;
    highp vec3 boxmax;
    int hit_next, miss_next;
};

group get_group(int which)
{
    group g;

    int j = which / group_data_row_size;
    int i = which - j * group_data_row_size;
    highp vec2 sample = vec2((float(i) + 0.25) / float(group_data_row_size), (float(j) + 0.25) / float(group_data_rows));

    g.boxmin = texture2D(group_boxmin, sample).xyz;
    g.boxmax = texture2D(group_boxmax, sample).xyz;

    highp vec2 group_child = texture2D(group_children, sample).xy;
    g.g1 = int(group_child.x);
    g.g2 = int(group_child.y);

    highp vec2 group_next = texture2D(group_hitmiss, sample).xy;
    g.hit_next = int(group_next.x);
    g.miss_next = int(group_next.y);

    highp vec2 group_object = texture2D(group_objects, sample).xy;
    g.start = int(group_object.x);
    g.count = int(group_object.y);

    g.is_branch = (g.count == 0);

    return g;
}

range group_bounds_intersect(in group g, in ray theray)
{
    return box_intersect(g.boxmin, g.boxmax, theray);
}

void sphere_intersect(int which, in ray theray, in range r, inout surface_hit hit)
{
    int j = which / sphere_data_row_size;
    int i = which - j * sphere_data_row_size;
    highp vec2 sample = vec2((float(i) + 0.25) / float(sphere_data_row_size), (float(j) + 0.25) / float(sphere_data_rows));

    highp vec4 sphere = texture2D(sphere_geometry, sample);
    highp vec3 center = sphere.xyz;
    float radius = sphere.w;

    float t = sphere_intersect(center, radius, theray, r);

    if(t == NO_t)
        return;

    if(t > hit.t)
        return;

    hit.which = which;
    hit.t = t;
}

void shade(in surface_hit hit, in ray theray, out highp vec4 normal, out highp vec4 point, out vec3 color)
{
    if(hit.which < 0) {
        color = hit.color;
        normal = vec4(0, 0, -1, 0);
        point = vec4(0, 0, 0, 1);
    } else {
        int j = hit.which / sphere_data_row_size;
        int i = hit.which - j * sphere_data_row_size;
        highp vec2 sample = vec2((float(i) + 0.25) / float(sphere_data_row_size), (float(j) + 0.25) / float(sphere_data_rows));

        highp vec4 sphere = texture2D(sphere_geometry, sample);
        highp vec4 center = vec4(sphere.x, sphere.y, sphere.z, 1.0);
        float radius = sphere.w;

        /*  snap to sphere surface */
        highp vec4 point0 = theray.o + theray.d * hit.t;
        highp vec4 to_surface = point0 - center;
        float distance = length(to_surface);

        point = center + to_surface * (radius / distance);
        normal = to_surface / radius;
        color = texture2D(sphere_colors, sample).xyz;
    }
}

/* limit to 200 for Windows Chrome < 30 */
const int max_bvh_iterations = 400;
const int max_leaf_tests = 8;

void group_intersect(int root, in ray theray, in range prevr, inout surface_hit hit)
{
    int g = root;

    /* 400 seems enough for the sorting pump */
    for(int i = 0; i < max_bvh_iterations; i++) {
        group gg = get_group(g);

        range r = range_intersect(prevr, group_bounds_intersect(gg, theray));

        if((!range_is_empty(r)) && (r.t0 < hit.t)) {
            // if(!gg.is_branch) {
                //  max 8 spheres in leaf, but limit on BVH depth
                //  takes precedence so could be fat leaves at max
                //  depth, need to carefully only make web scenes with
                //  8 or fewer spheres at leaf
                for(int j = 0; j < max_leaf_tests; j++) {
                    if(j >= gg.count)
			break;
		    sphere_intersect(gg.start + j, theray, r, hit);
		}
            // }
            g = gg.hit_next;
        } else {
            g = gg.miss_next;
        }

        if(g >= 0xffff)
            return;
        if(i == max_bvh_iterations - 1)
            set_bad_hit(hit, 1.0, 1.0, 0.0);
    }
}

highp vec3 eye_ray(float u, float v, float aspect, float fov)
{
    float eye_alpha = fov * (u - 0.5);
    float eye_beta = fov * (v - 0.5) * aspect;

    /* eye space ray */
    highp vec3 eye;
    eye.x = sin(eye_alpha) * cos(eye_beta);
    eye.y = sin(eye_beta);
    eye.z = cos(eye_alpha) * cos(eye_beta);

    return eye;
}

const bool cast_shadows = true;

void transform(in ray r, in mat4 matrix, in mat4 normal_matrix, out ray t)
{
    t.o = matrix * r.o;
    t.d = normal_matrix * r.d;
}

bool do_one_whitted(in ray worldray, out vec3 result, out ray reflected)
{
    surface_hit shading = surface_hit_init(background_color);

    ray objectray;
    transform(worldray, object_matrix, object_normal_matrix, objectray);

    group_intersect(tree_root, objectray, range(0.0, 100000000.0), shading);

    if(shading.t >= infinitely_far) {

        result = background_color;
        return false;

    } else {

        vec3 color;
        highp vec4 object_normal, object_point;
        highp vec4 world_normal, world_point;
        shade(shading, objectray, object_normal, object_point, color);
        world_normal = object_normal_inverse * object_normal;
        world_point = object_inverse * object_point;

	ray world_shadowray;
	ray object_shadowray;

	world_shadowray.o = world_point + world_normal * 0.0001;
	world_shadowray.d = vec4(light_dir.x, light_dir.y, light_dir.z, 0);
	surface_hit shadow_hit = surface_hit_init(vec3(0.0, 0.0, 0.0));
        transform(world_shadowray, object_matrix, object_normal_matrix, object_shadowray);
	group_intersect(tree_root, object_shadowray, range(0.0, 100000000.0), shadow_hit);

        if(shadow_hit.t < infinitely_far) {
            result = color * 0.1;
        } else {
            float diffuse = max(0.1, dot(world_normal, light_dir));
            result = color * diffuse;
        }

        reflected.o = world_point + world_normal * .0001;
        reflected.d = reflect(worldray.d, world_normal);
        return true;
    }
}

const int bounce_count = 2;

void main()
{
    gl_FragColor = vec4(0, 0, 0, 1);

    ray eyeray, worldray, objectray;

    eyeray.d.xyz = eye_ray(ftex.x, 1.0 - ftex.y, aspect, fov);
    eyeray.d.w = 0.0;
    eyeray.o = vec4(0.0, 0.0, 0.0, 1.0);

    transform(eyeray, camera_matrix, camera_normal_matrix, worldray);

    float intensity = 1.0;
    for(int i = 0; i < bounce_count; i++) {
        ray reflected;
        vec3 color;
        bool hit_something = do_one_whitted(worldray, color, reflected);
        gl_FragColor.xyz += color * intensity;
        if(!hit_something)
            break;
        intensity *= .5;
        worldray = reflected;
    }
}
