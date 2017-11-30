#include "tiny_obj_loader.h"
#include "ppm_example.c"
#include "vec4.h"
#include "mat4.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// triangle struct
struct face {
    vec4 vertex1;
    vec4 vertex2;
    vec4 vertex3;
};


float l, r, t, b, n, f, eye_x, eye_y, eye_z, center_x, center_y, center_z, up_x, up_y, up_z;

vector<tinyobj::shape_t> shapes;
vector<tinyobj::material_t> materials;


int width;
int height;
char* outputPpm;
char* colorOption;
bool isColorOption;

void readCam(char* camText);
mat4 orientMat();
void rasterizeTriangles();
void colorPixel(int pixel, int numTri);


img_t *img;
vector<vec4> normVecs;
vector<face> faceVecs;
vector<float> indices;
vector<vec4> transNormVecs;

int main(int argc, char *argv[]) {

    // inputs
    char* inputFile = argv[1];
    char* cameraText = argv[2];
    width = atoi(argv[3]);
    height = atoi(argv[4]);
    outputPpm = argv[5];

    if (argc == 7) {
        colorOption = argv[6];
        isColorOption = true;
    }

    std::string path = "/Users/catherineyang/Documents/Sophomore Spring 2017/CIS 460/new_rasterizer/rasterizer\\";
    const char *path1 = path.c_str();
    cout<< tinyobj::LoadObj(shapes, materials, inputFile);

    readCam(cameraText);

    img = new_img(width, height);

    rasterizeTriangles();

    return 0;
}

// read from camera file
void readCam(char *camText) {
    FILE *camFile = fopen(camText, "r");
    fscanf(camFile, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
           &l, &r, &t, &b, &n, &f, &eye_x, &eye_y, &eye_z, &center_x, &center_y, &center_z, &up_x, &up_y, &up_z);
    fclose(camFile);

}

// orientatio matrix
mat4 orientMat() {
    vec4 up = vec4(up_x, up_y, up_z, 0).normalize();

    vec4 f = vec4(center_x - eye_x, center_y - eye_y, center_z - eye_z, 0).normalize();
    vec4 r = cross(f, up).normalize();
    vec4 u = cross(f, r).normalize();

    mat4 orientation = mat4(vec4(r[0], u[0], f[0], 0),
                            vec4(r[1], u[1], f[1], 0),
                            vec4(r[2], u[2], f[2], 0),
                            vec4(0, 0, 0, 1));
    return orientation;
}

// view matrix
mat4 viewMat() {
    mat4 translation = mat4(vec4(1, 0, 0, 0),
                            vec4(0, 1, 0, 0),
                            vec4(0, 0, 1, 0),
                            vec4(-eye_x, -eye_y, -eye_z, 1));
    return orientMat() * translation;
}

// projection matrix
mat4 projMat() {
    mat4 frustum = mat4(vec4((2*n)/(r-l), 0, 0, 0),
                        vec4(0, (2*n)/(t-b), 0, 0),
                        vec4((r+l)/(r-l), (t+b)/(t-b), f/(f-n), 1),
                        vec4(0, 0, (-f*n)/(f-n), 0));

    return frustum;
}

void rasterizeTriangles() {

    float zBuffs[width * height];

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            zBuffs[width * i + j] = 2;
        }
    }

    // positions: vector of vertex coordinates
    vector<float> positions;
    for (int i = 0; i < static_cast<int>(shapes.size()); i++) {
        for (int j = 0; j < static_cast<int>(shapes[i].mesh.positions.size()); j++) {
            positions.push_back(shapes[i].mesh.positions[j]);
        }
    }

    // vector of vec4s, with w coordinate of 1
    vector<vec4> posVecs;
    for (int i = 0; i < static_cast<int>(positions.size()); i+=3) {
        vec4 here = vec4(positions[i], positions[i+1], positions[i+2], 1);
        posVecs.push_back(here);
    }

    // normals; vector of normal coordinates
    vector<float> normals;

    for (int i = 0; i < static_cast<int>(shapes.size()); i++) {
        for (int j = 0; j < static_cast<int>(shapes[i].mesh.normals.size()); j++) {
            normals.push_back(shapes[i].mesh.normals[j]);
        }
    }

    // vector of vec4s, with w coordinate of 0
    for (int i = 0; i < static_cast<int>(normals.size()); i+=3) {
        vec4 here = vec4(normals[i], normals[i+1], normals[i+2], 0);
        normVecs.push_back(here);
    }

    for (int i = 0; i < static_cast<int>(normVecs.size()); i++) {
        vec4 update = orientMat() * normVecs[i];
        transNormVecs.push_back(update);
    }

    // indices: vector of vertex indices that specify vertices of triangles
    for (int i = 0; i < static_cast<int>(shapes.size()); i++) {
        for (int j = 0; j < static_cast<int>(shapes[i].mesh.indices.size()); j++) {
            indices.push_back(shapes[i].mesh.indices[j]);
        }
    }

    // vector of faces
    for (int i = 0; i < static_cast<int>(indices.size()); i+=3) {
        face here;
        here.vertex1 = posVecs[indices[i]];
        here.vertex2 = posVecs[indices[i+1]];
        here.vertex3 = posVecs[indices[i+2]];
        faceVecs.push_back(here);
    }

    // calculations
    for (int i = 0; i < static_cast<int>(faceVecs.size()); i++) {
        vec4 first = projMat() * viewMat() * faceVecs[i].vertex1;
        float w1 = first[3];
        faceVecs[i].vertex1 = first/w1;

        vec4 second = projMat() * viewMat() * faceVecs[i].vertex2;
        float w2 = second[3];
        faceVecs[i].vertex2 = second/w2;

        vec4 third = projMat() * viewMat() * faceVecs[i].vertex3;
        float w3 = third[3];
        faceVecs[i].vertex3 = third/w3;
    }

    for (int i = 0; i < static_cast<int>(faceVecs.size()); i++) {
        // convert from NDC to pixel coordinates
        vec4 first = faceVecs[i].vertex1;
        float x1 = (first[0] + 1) * (width/2);
        float y1 = (1 - first[1]) * (height/2);
        float z1 = first[2];
        float w1 = first[3];
        faceVecs[i].vertex1 = vec4(x1, y1, z1, w1);

        vec4 second = faceVecs[i].vertex2;
        float x2 = (second[0] + 1) * (width/2);
        float y2 = (1 - second[1]) * (height/2);
        float z2 = second[2];
        float w2 = second[3];
        faceVecs[i].vertex2 = vec4(x2, y2, z2, w2);

        vec4 third = faceVecs[i].vertex3;
        float x3 = (third[0] + 1) * (width/2);
        float y3 = (1 - third[1]) * (height/2);
        float z3 = third[2];
        float w3 = third[3];
        faceVecs[i].vertex3 = vec4(x3, y3, z3, w3);


        // triangles that are not visible
        if (z1 < 0 && z2 < 0 && z3 < 0) {
            continue;
        } else if (z1 > 1 && z2 > 1 && z2 > 1) {
            continue;
        }

        // compute 2D bounding box
        int minX = floor(std::min({x1, x2, x3}));
        int maxX = ceil(std::max({x1, x2, x3}));

        int minY = floor(std::min({y1, y2, y3}));
        int maxY = ceil(std::max({y1, y2, y3}));

        if (maxX < 0 || minX > width || maxY < 0 || minY > height) {
            continue;
        }

        // vertex 1 and 2
        float m12 = (faceVecs[i].vertex2[1] - faceVecs[i].vertex1[1]) /
                    (faceVecs[i].vertex2[0] - faceVecs[i].vertex1[0]);

        // vertex 1 and 3
        float m13 = (faceVecs[i].vertex3[1] - faceVecs[i].vertex1[1]) /
                    (faceVecs[i].vertex3[0] - faceVecs[i].vertex1[0]);

        // vertex 2 and 3
        float m23 = (faceVecs[i].vertex3[1] - faceVecs[i].vertex2[1]) /
                    (faceVecs[i].vertex3[0] - faceVecs[i].vertex2[0]);


        float enter12, enter13, enter23;
        float exit12, exit13, exit23;

        for (float y = (minY + 0.5); y <= maxY; y++) {
            // x values that are outside appropriate range
            enter12 = width * 2;
            enter13 = width * 2;
            enter23 = width * 2;
            exit12 = -1;
            exit13 = -1;
            exit23 = -1;

            // checking for intersections for each of the 3 sides of the triangle
            if (m12 == 0) { // horizontal
                if (faceVecs[i].vertex2[1] <= y + 0.5 && faceVecs[i].vertex2[1] >= y - 0.5) {
                    enter12 = std::min(faceVecs[i].vertex2[0], faceVecs[i].vertex1[0]);
                    exit12 = std::max(faceVecs[i].vertex2[0], faceVecs[i].vertex1[0]);
                }
            } else if (m12 == INFINITY || m12 == -INFINITY) {
                float startY = std::min(faceVecs[i].vertex1[1], faceVecs[i].vertex2[1]);
                float endY = std::max(faceVecs[i].vertex1[1], faceVecs[i].vertex2[1]);
                if (y <= endY && y >= startY) {
                    enter12 = faceVecs[i].vertex2[0];
                    exit12 = faceVecs[i].vertex2[0];
                }
            } else {
                float x12 = (y + (m12 * faceVecs[i].vertex2[0]) - faceVecs[i].vertex2[1]) / m12;
                float startX = std::min(faceVecs[i].vertex2[0], faceVecs[i].vertex1[0]);
                float endX = std::max(faceVecs[i].vertex2[0], faceVecs[i].vertex1[0]);
                if (x12 <= endX && x12 >= startX) {
                    enter12 = x12;
                    exit12 = x12;
                }
            }


            if (m13 == 0) { // horizontal
                if (faceVecs[i].vertex3[1] <= y + 0.5 && faceVecs[i].vertex3[1] >= y - 0.5) {
                    enter13 = std::min(faceVecs[i].vertex3[0], faceVecs[i].vertex1[0]);
                    exit13 = std::max(faceVecs[i].vertex3[0], faceVecs[i].vertex1[0]);
                }
            } else if (m13 == INFINITY || m13 == -INFINITY) {
                float startY = std::min(faceVecs[i].vertex3[1], faceVecs[i].vertex1[1]);
                float endY = std::max(faceVecs[i].vertex3[1], faceVecs[i].vertex1[1]);
                if (y <= endY && y >= startY) {
                    enter13 = faceVecs[i].vertex3[0];
                    exit13 = faceVecs[i].vertex3[0];
                }
            } else {
                float x13 = (y + (m13 * faceVecs[i].vertex3[0]) - faceVecs[i].vertex3[1]) / m13;
                float startX = std::min(faceVecs[i].vertex3[0], faceVecs[i].vertex1[0]);
                float endX = std::max(faceVecs[i].vertex3[0], faceVecs[i].vertex1[0]);
                if (x13 <= endX && x13 >= startX) {
                    enter13 = x13;
                    exit13 = x13;
                }
            }


            if (m23 == 0) { // horizontal
                if (faceVecs[i].vertex3[1] <= y + 0.5 && faceVecs[i].vertex3[1] >= y - 0.5) {
                    enter23 = std::min(faceVecs[i].vertex3[0], faceVecs[i].vertex2[0]);
                    exit23 = std::max(faceVecs[i].vertex3[0], faceVecs[i].vertex2[0]);
                }
            } else if (m23 == INFINITY || m23 == -INFINITY) {
                float startY = std::min(faceVecs[i].vertex3[1], faceVecs[i].vertex2[1]);
                float endY = std::max(faceVecs[i].vertex3[1], faceVecs[i].vertex2[1]);
                if (y <= endY && y >= startY) {
                    enter23 = faceVecs[i].vertex3[0];
                    exit23 = faceVecs[i].vertex3[0];
                }
            } else {
                float x23 = (y + (m23 * faceVecs[i].vertex3[0]) - faceVecs[i].vertex3[1]) / m23;
                float startX = std::min(faceVecs[i].vertex3[0], faceVecs[i].vertex2[0]);
                float endX = std::max(faceVecs[i].vertex3[0], faceVecs[i].vertex2[0]);
                if (x23 <= endX && x23 >= startX) {
                   enter23 = x23;
                   exit23 = x23;
                }
            }

            int enter = floor(std::min(enter12, std::min(enter13, enter23)));
            int exit = ceil(std::max(exit12, std::max(exit13, exit23)));

            float x1 = faceVecs[i].vertex1[0];
            float y1 = faceVecs[i].vertex1[1];
            float x2 = faceVecs[i].vertex2[0];
            float y2 = faceVecs[i].vertex2[1];
            float x3 = faceVecs[i].vertex3[0];
            float y3 = faceVecs[i].vertex3[1];

            for (int x = enter; x < exit; x++) {

                float s = 0.5 * cross(vec4(x1, y1, 0, 0) - vec4(x2, y2, 0, 0), vec4(x3, y3, 0, 0) - vec4(x2, y2, 0, 0)).length();
                float s1 = 0.5 * cross(vec4(x, y, 0, 0) - vec4(x2, y2, 0, 0), vec4(x, y, 0, 0) - vec4(x3, y3, 0, 0)).length();
                float s2 = 0.5 * cross(vec4(x, y, 0, 0) - vec4(x3, y3, 0, 0), vec4(x, y, 0, 0) - vec4(x1, y1, 0, 0)).length();
                float s3 = 0.5 * cross(vec4(x, y, 0, 0) - vec4(x1, y1, 0, 0), vec4(x, y, 0, 0) - vec4(x2, y2, 0, 0)).length();

                float w1 = s1 / s;
                float w2 = s2 / s;
                float w3 = s3 / s;
                float Z, pz1, pz2;


                int pixel = width * (int)y + x;

                float alpha1, alpha2;

                if ((y2 - y1) == 0) {
                    alpha1 = (x - x1) / (x2 - x1);
                    alpha2 = (y - y1) / (y3 - y1);
                } else if ((y3 - y1) == 0) {
                    alpha1 = (y - y1) / (y2 - y1);
                    alpha2 = (x - x1) / (x3 - x1);
                } else {
                    alpha1 = (y - y1) / (y2 - y1);
                    alpha2 = (y - y1) / (y3 - y1);
                }

                float alpha3 = (x - enter) / (exit - enter);

                if (strcmp(colorOption, "--norm_gouraud_z") == 0) {
                    pz1 = 1 / (((1 - alpha1) / faceVecs[i].vertex1[2]) + (alpha1 / faceVecs[i].vertex2[2]));
                    pz2 = 1 / (((1 - alpha2) / faceVecs[i].vertex1[2]) + (alpha2 / faceVecs[i].vertex3[2]));
                    Z = 1 / (((1 - alpha3) / pz1) + (alpha3 / pz2));
                } else if (strcmp(colorOption, "--norm_bary_z") == 0) {
                    w1 = w1 / faceVecs[i].vertex1[2];
                    w2 = w2 / faceVecs[i].vertex2[2];
                    w3 = w3 / faceVecs[i].vertex3[2];
                    Z = 1 / (w1 + w2 + w3);
                } else {
                    w1 = w1 * faceVecs[i].vertex1[2];
                    w2 = w2 * faceVecs[i].vertex2[2];
                    w3 = w3 * faceVecs[i].vertex3[2];
                    Z = w1 + w2 + w3;
                }

                vec4 norm1 = transNormVecs[indices[i * 3]];
                vec4 norm2 = transNormVecs[indices[(i * 3) + 1]];
                vec4 norm3 = transNormVecs[indices[(i * 3) + 2]];

                float R1, G1, B1;

                // normalized gouraud shading
                if (strcmp(colorOption, "--norm_gouraud") == 0) {
                    float pR1 = ((1.0f - alpha1) * norm1[0]) + (alpha1 * norm2[0]);
                    float pR2 = ((1.0f - alpha2) * norm1[0]) + (alpha2 * norm3[0]);
                    float pG1 = ((1.0f - alpha1) * norm1[1]) + (alpha1 * norm2[1]);
                    float pG2 = ((1.0f - alpha2) * norm1[1]) + (alpha2 * norm3[1]);
                    float pB1 = ((1.0f - alpha1) * norm1[2]) + (alpha1 * norm2[2]);
                    float pB2 = ((1.0f - alpha2) * norm1[2]) + (alpha2 * norm3[2]);

                    R1 = ((1.0f - alpha3) * pR1) + (alpha3 * pR2);
                    G1 = ((1.0f - alpha3) * pG1) + (alpha3 * pG2);
                    B1 = ((1.0f - alpha3) * pB1) + (alpha3 * pB2);
                } else if (strcmp(colorOption, "--norm_gouraud_z") == 0) { 
                    float pR1 = 1 / (((1 - alpha1) / norm1[0]) + (alpha1 / norm2[0]));
                    float pR2 = 1 / (((1 - alpha2) / norm1[0]) + (alpha2 / norm3[0]));
                    float pG1 = 1 / (((1 - alpha1) / norm1[1]) + (alpha1 / norm2[1]));
                    float pG2 = 1 / (((1 - alpha2) / norm1[1]) + (alpha2 / norm3[1]));
                    float pB1 = 1 / (((1 - alpha1) / norm1[2]) + (alpha1 / norm2[2]));
                    float pB2 = 1 / (((1 - alpha2) / norm1[2]) + (alpha2 / norm3[2]));
                    // problematic because some z values of the normals are 0 (dividing by 0)

                    R1 = 1 / (((1 - alpha3) / pR1) + (alpha3 / pR2));
                    G1 = 1 / (((1 - alpha3) / pG1) + (alpha3 / pG2));
                    B1 = 1 / (((1 - alpha3) / pB1) + (alpha3 / pB2));
                } else if (strcmp(colorOption, "--norm_bary") == 0) { // barycentric interpolation
                    R1 = w1 * norm1[0] + w2 * norm2[0] + w3 * norm3[0];
                    G1 = w1 * norm1[1] + w2 * norm2[1] + w3 * norm3[1];
                    B1 = w1 * norm1[2] + w2 * norm2[2] + w3 * norm3[2];
                } else if (strcmp(colorOption, "--norm_bary_z") == 0) {
                    R1 = Z * (w1 * norm1[0] + w2 * norm2[0] + w3 * norm3[0]);
                    G1 = Z * (w1 * norm1[1] + w2 * norm2[1] + w3 * norm3[1]);
                    B1 = Z * (w1 * norm1[2] + w2 * norm2[2] + w3 * norm3[2]);
                }

                if (Z > 0 && Z < 1 && Z < zBuffs[pixel]) {
                    if ((strcmp(colorOption, "--norm_gouraud") == 0) || (strcmp(colorOption, "--norm_bary") == 0) ||
                        (strcmp(colorOption, "--norm_bary_z") == 0) || (strcmp(colorOption, "--norm_gouraud_z") == 0)) {
                        float R = (255 / 2.0f) * (R1 + 1.0f);
                        float G = (255 / 2.0f) * (G1 + 1.0f);
                        float B = (255 / 2.0f) * (B1 + 1.0f);
                        float r = (int)R;
                        float g = (int)G;
                        float b = (int)B;

                        if (r < 0) { r = 0;}
                        if (g < 0) { g = 0;}
                        if (b < 0) { b = 0;}
                        if (r > 255) { r = 255;}
                        if (g > 255) { g = 255;}
                        if (b > 255) { b = 255;}

                        img->data[pixel].r = r;
                        img->data[pixel].g = g;
                        img->data[pixel].b = b;
                    } else {
                        colorPixel(pixel, i);
                    }
                    zBuffs[pixel] = Z;
                }

            }
        }
    }

    write_ppm(img, outputPpm);
}


void colorPixel(int pixel, int numTri) {
    int r, g, b;

    if (!isColorOption) {
        float color[3];
        for (int i = 0; i < static_cast<int>(materials.size()); i++) {
            for (int j = 0; j < 3; j++) {
                color[j] = materials[i].diffuse[j];
            }
        }
        r = color[0] * 255;
        g = color[1] * 255;
        b = color[2] * 255;
    } else if (strcmp(colorOption, "--white") == 0) {
        r = 255;
        g = 255;
        b = 255;
    } else if (strcmp(colorOption, "--norm_flat") == 0) {
        float R = (255 / 2.0f) * (transNormVecs[indices[numTri * 3]][0] + 1.0f);
        float G = (255 / 2.0f) * (transNormVecs[indices[numTri * 3]][1] + 1.0f);
        float B = (255 / 2.0f) * (transNormVecs[indices[numTri * 3]][2] + 1.0f);

        r = (int)R;
        g = (int)G;
        b = (int)B;
    }

    img->data[pixel].r = r;
    img->data[pixel].g = g;
    img->data[pixel].b = b;
}

