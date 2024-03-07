#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <cglm/cglm.h>

typedef struct ray ray;
typedef struct hit_record hit_record;
typedef struct sphere sphere;
typedef bool (*hit_fp)(ray *r, float tmin, float tmax, void *obj,
                       hit_record *hr);
typedef struct hittable_list hittable_list;
typedef struct camera camera;

struct ray
{
    vec3 origin;
    vec3 direction;
};

struct hit_record
{
    vec3 position;
    vec3 normal;
    float t;
    bool front_face;
};

struct sphere
{
    vec3 center;
    float radius;
};

struct hittable_list
{
    hit_fp *hit_funcs;
    void **hittables;
    int length;
    int capacity;
};

struct camera
{
    int image_width;
    int image_height;
    float aspect_ratio;
    vec3 center;
    vec3 pixel00_loc;
    vec3 pixel_delta_u;
    vec3 pixel_delta_v;
};

/**
 * Calculates the point of intersection (if any) of the ray with a hittable
 *
 * @param[in] r The ray
 * @param[in] t The length along the ray of the intersection
 * @param[out] dest The vec3 location of the intersection in world space
 */
void ray_at(ray *r, float t, float *dest);

/**
 * Determines the color of the pixel
 *
 * @param[in] r The ray used to find the pixel color
 * @param[out] dest The resultant pixel color
 */
void ray_color(ray *r, hittable_list *hl, float *dest);

/**
 * Determines if the ray intersects with the sphere
 *
 * @param[in] r The ray
 * @param[in] tmin The minimum value along the ray to calculate intersection
 * @param[in] tmax The maximum value along the ray to calculate intersection
 * @param[in] obj The type of object with which the ray intersects
 * @param[out] hr The intersection information
 *
 * @return Whether or not the ray intersects with the object
 */
bool hit_sphere(ray *r, float tmin, float tmax, void *obj, hit_record *hr);

/**
 * Chooses which object normal should be used for the hit record
 *
 * @param[in] r The ray
 * @param[in] outward_normal The outward facing normal vector
 * @param[out] hr The intersection information
 */
void set_face_normal(ray *r, float *outward_normal, hit_record *hr);

/**
 * Initializes a list of hittable objects within the scene
 *
 * @return The initialized hittable list
 */
hittable_list *create_hittable_list(void);

/**
 * Deletes a list of hittable objects
 *
 * @param[in] hl The list of hittable objects
 */
void delete_hittable_list(hittable_list *hl);

/**
 * Adds a hittable object to the end of the list
 *
 * @param[in] hl The list of hittable objects
 * @param[in] hfp The ray-object intersection function pointer
 * @param[in] obj The object
 */
void hittable_list_add(hittable_list *hl, hit_fp hfp, void *obj);

/**
 * Removes a hittable object from the end of the list
 * 
 * @param[in] hl The list of hittable objects
 */
void hittable_list_remove(hittable_list *hl);

/**
 * Initializes the scene camera
 *
 * @return The intialized camera
 */
camera *create_camera(void);

/**
 * Deletes the scene camera
 *
 * @param[in] cam The camera to delete
 */
void delete_camera(camera *cam);

/**
 * Renders the scene
 *
 * @param[in] hl The list of hittable objects
 * @param[in] camera The scene camera
 */
void render(hittable_list *hl, camera *cam);

int
main(void)
{
    hittable_list *world = create_hittable_list();
    camera *cam = create_camera();

    hittable_list_add(world, hit_sphere,
                      (void *)&(sphere){{0.0f, 0.0f, -1.0f}, 0.5f});
    hittable_list_add(world, hit_sphere,
                      (void *)&(sphere){{0.0f, -100.5f, -1.0f}, 100.0f});

    render(world, cam);
    printf("Done\n");

    free(cam);
    delete_hittable_list(world);

    return 0;
}

void
ray_at(ray *r, float t, float *dest)
{
    vec3 temp;

    glm_vec3_scale(r->direction, t, temp);
    glm_vec3_add(r->origin, temp, dest);
}

void
ray_color(ray *r, hittable_list *hl, float *dest)
{
    int i;
    float a;
    vec3 temp;
    vec3 unit_direction;
    hit_record hr = {GLM_VEC3_ZERO_INIT, GLM_VEC3_ZERO_INIT, 0.0f, true};

    for (i = 0; i < hl->length; i++) {
        if (hl->hit_funcs[i](r, 0.0f, INFINITY, hl->hittables[i], &hr)) {
            glm_vec3_add(hr.normal, GLM_VEC3_ONE, temp);
            glm_vec3_scale(temp, 0.5f, dest);
            return;
        }
    }

    glm_vec3_normalize_to(r->direction, unit_direction);
    a = 0.5f * (unit_direction[1] + 1.0f);
    glm_vec3_scale(GLM_VEC3_ONE, 1.0f - a, dest);
    glm_vec3_muladds((vec3){0.5f, 0.7f, 1.0f}, a, dest);
}

bool
hit_sphere(ray *r, float tmin, float tmax, void *obj, hit_record *hr)
{
    float discriminant;
    float sqd, root;
    float a, b, c;
    vec3 oc, temp, outward_normal;
    sphere *s = (sphere *)obj;

    glm_vec3_sub(r->origin, s->center, oc);

    a = glm_vec3_dot(r->direction, r->direction);
    b = glm_vec3_dot(oc, r->direction);
    c = glm_vec3_dot(oc, oc) - s->radius * s->radius;

    discriminant = b * b - a * c;
    if (discriminant < 0)
        return false;
    
    sqd = sqrtf(discriminant);
    root = (-b - sqd) / a;

    if (root <= tmin || root >= tmax) {
        root = (-b + sqd) / a;
        if (root < tmin || root >= tmax)
            return false;
    }

    hr->t = root;
    ray_at(r, hr->t, hr->position);
    glm_vec3_sub(hr->position, s->center, temp);
    glm_vec3_divs(temp, s->radius, outward_normal);
    set_face_normal(r, outward_normal, hr);

    return true;
}

void
set_face_normal(ray *r, float *outward_normal, hit_record *hr)
{
    vec3 outward_normal_neg;
    glm_vec3_negate_to(outward_normal, outward_normal_neg);

    hr->front_face = glm_vec3_dot(r->direction, outward_normal) < 0.0f;

    hr->front_face ? glm_vec3_copy(outward_normal, hr->normal)
                   : glm_vec3_copy(outward_normal_neg, hr->normal);
}

hittable_list *
create_hittable_list(void)
{
    hittable_list *hl = malloc(sizeof(*hl));

    if (hl == NULL) {
        fprintf(stderr, "Error: Could not allocate enough memory for hittable "
                "list\n");
        exit(EXIT_FAILURE);
    }

    hl->capacity = 8;
    hl->length = 0;

    hl->hit_funcs = calloc(hl->capacity, sizeof(*hl->hit_funcs));
    if (hl->hit_funcs == NULL) {
        fprintf(stderr, "Error: Could not allocate enough memory for list of "
                "function pointers in hittable list\n");
        exit(EXIT_FAILURE);
    }

    hl->hittables = calloc(hl->capacity, sizeof(*hl->hittables));
    if (hl->hittables == NULL) {
        fprintf(stderr, "Error: Could not allocate enough memory for list of "
                "hittables in hittable list\n");
        exit(EXIT_FAILURE);
    }

    return hl;
}

void
delete_hittable_list(hittable_list *hl)
{
    free(hl->hittables);
    free(hl->hit_funcs);
    free(hl);
}

void
hittable_list_add(hittable_list *hl, hit_fp hfp, void *obj)
{
    hit_fp *temp_hfp;
    void **temp_vp;

    if (hl->length - 1 == hl->capacity) {
        temp_hfp = realloc(hl->hit_funcs,
                           hl->capacity * 2 * sizeof(*hl->hit_funcs));

        if (temp_hfp == NULL) {
            fprintf(stderr, "Error: Failed to increase space for hittable list "
                    "function pointer array\n");
            delete_hittable_list(hl);
            exit(EXIT_FAILURE);
        } else {
            hl->hit_funcs = temp_hfp;
        }

        temp_vp = realloc(hl->hittables,
                          hl->capacity * 2 * sizeof(*hl->hittables));

        if (temp_vp == NULL) {
            fprintf(stderr, "Error: Failed to increase space for hittable list "
                    "hittables array\n");
            delete_hittable_list(hl);
            exit(EXIT_FAILURE);
        } else {
            hl->hittables = temp_vp;
        }

        hl->capacity *= 2;
    }

    hl->hit_funcs[hl->length] = hfp;
    hl->hittables[hl->length] = obj;
    hl->length++;
}

void
hittable_list_remove(hittable_list *hl)
{
    hit_fp *temp_hfp;
    void **temp_vp;

    hl->hit_funcs[hl->length] = NULL;
    hl->hittables[hl->length] = NULL;

    if (hl->length > 0)
        hl->length--;

    if (hl->length == hl->capacity / 4 && hl->capacity > 8) {
        temp_hfp = realloc(hl->hit_funcs,
                           hl->capacity / 2 * sizeof(*hl->hit_funcs));

        if (temp_hfp == NULL) {
            fprintf(stderr, "Error: Failed to decrease space for hittable list "
                    "function pointer array\n");
            delete_hittable_list(hl);
            exit(EXIT_FAILURE);
        } else {
            hl->hit_funcs = temp_hfp;
        }

        temp_vp = realloc(hl->hittables,
                          hl->capacity / 2 * sizeof(*hl->hittables));

        if (temp_vp == NULL) {
            fprintf(stderr, "Error: Failed to decrease space for hittable list "
                    "hittables array\n");
            delete_hittable_list(hl);
            exit(EXIT_FAILURE);
        } else {
            hl->hittables = temp_vp;
        }

        hl->capacity /= 2;
    }
}

camera *
create_camera(void)
{
    float focal_length = 1.0f;
    float viewport_width, viewport_height;
    vec3 viewport_u = GLM_VEC3_ZERO_INIT;
    vec3 viewport_v = GLM_VEC3_ZERO_INIT;
    vec3 viewport_upper_left, temp;

    camera *cam = malloc(sizeof(*cam));

    if (cam == NULL) {
        fprintf(stderr, "Error: Could not allocate enough memory for camera\n");
        exit(EXIT_FAILURE);
    }

    glm_vec3_zero(cam->center);

    /* ----- Image ----- */
    cam->aspect_ratio = 16.0f / 9.0f;
    cam->image_width = 400;

    cam->image_height = (int)(cam->image_width / cam->aspect_ratio);
    cam->image_height = (cam->image_height < 1) ? 1 : cam->image_height;

    /* ----- Camera ----- */
    viewport_height = 2.0f;
    viewport_width = viewport_height * ((float)cam->image_width / cam->image_height);

    viewport_u[0] = viewport_width;
    viewport_v[1] = -viewport_height;

    glm_vec3_divs(viewport_u, cam->image_width, cam->pixel_delta_u);
    glm_vec3_divs(viewport_v, cam->image_height, cam->pixel_delta_v);

    /* I wish there were a divsubs function like mulsubs */
    glm_vec3_sub(cam->center, (vec3){0, 0, focal_length},
                 viewport_upper_left);
    glm_vec3_mulsubs(viewport_u, 0.5f, viewport_upper_left);
    glm_vec3_mulsubs(viewport_v, 0.5f, viewport_upper_left);

    glm_vec3_add(cam->pixel_delta_u, cam->pixel_delta_v, temp);
    glm_vec3_scale(temp, 0.5f, temp);
    glm_vec3_add(viewport_upper_left, temp, cam->pixel00_loc);

    return cam;
}

void
delete_camera(camera *cam)
{
    free(cam);
}

void
render(hittable_list *hl, camera *cam)
{
    int i, j;
    int ir, ig, ib;
    vec3 pixel_color, pixel_center, temp;

    ray r = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

    FILE *img_fp = NULL;
    const char *filename = "image.ppm";

    if ((img_fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(img_fp, "P3\n%d %d\n255\n", cam->image_width, cam->image_height);

    for (i = 0; i < cam->image_height; i++) {
        for (j = 0; j < cam->image_width; j++) {
            glm_vec3_scale(cam->pixel_delta_u, j, temp);
            glm_vec3_add(cam->pixel00_loc, temp, pixel_center);
            glm_vec3_muladds(cam->pixel_delta_v, i, pixel_center);

            glm_vec3_sub(pixel_center, cam->center, r.direction);

            ray_color(&r, hl, pixel_color);

            ir = (int)(255.99f * pixel_color[0]);
            ig = (int)(255.99f * pixel_color[1]);
            ib = (int)(255.99f * pixel_color[2]);

            fprintf(img_fp, "%d %d %d\n", ir, ig, ib);
        }
    }

    fclose(img_fp);
}
/* EOF */
