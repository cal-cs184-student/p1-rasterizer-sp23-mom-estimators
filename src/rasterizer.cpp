#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

        sample_buffer[y * width + x] = c;
//    rgb_framebuffer_target[3 * (y * width + x)] = (unsigned char)(c.r * 255);
//    rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(c.g * 255);
//    rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(c.b * 255);

  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.

  // CREATED FOR TASK 1: indicator function
  // check if a point is inside the triangle
  // maggie's shitty code that doesn't make sense
//   this is a little buggy but almost works?
  bool is_point_in_triangle(float x0 , float y0, float x1, float y1, float x2, float y2, float ptx, float pty) {
      // three line test
      float pts[3][2] = {{x0, y0}, {x1, y1} , {x2, y2}};
      float dx[] = {abs(x1-x0), abs(x2-x0), abs(x2-x1)};
      float dy[] = {abs(y1-y0), abs(y2-y0), abs(y2-y1)};
      float l1[] = {0, 0, 0};
      float l2[] = {0, 0, 0};
      float l3[] = {0, 0, 0};

//      std::cout << "hi" << std::endl;

      for (int i = 0; i < 3; i++) {
          if (i == 0) {
              l1[0] = -1 * dy[0]; // A_i x
              l1[1] = dx[0]; // B_i y
              l1[2] = pts[i][0] * dy[i] - pts[i][1] * dx[i]; // C_i
          };
          if (i == 1) {
              l2[0] = -1 * dy[0]; // A_i
              l2[1] = dx[0]; // B_i
              l2[2] = pts[i][0] * dy[i] - pts[i][1] * dx[i]; // C_i
          };
          if (i == 2) {
              l3[0] = -1 * dy[0]; // A_i
              l3[1] = dx[0]; // B_i
              l3[2] = pts[i][0] * dy[i] - pts[i][1] * dx[i]; // C_i
          };
      }

      if ((l1[0] * ptx + l1[1] * pty + l1[2]) >= 0 && (l2[0] * ptx + l2[1] * pty + l2[2]) >= 0 && (l3[0] * ptx + l3[1] * pty + l3[2]) >= 0) {
          return 1;
      }
      return 0;
  }

//  using the cross product / normal vectors to determine if it is inside the triangle
  bool is_in_triangle(float x0, float x1, float y0, float y1, float x2, float y2, float x, float y) {
      float x_coords[] = {x0, x1, x2};
      float y_coords[] = {y0, y1, y2};
//      float crosses[3];
//
//      for (int i = 0; i < 3; i++) {
//          Vector2D edge_vec = Vector2D(x_coords[i], y_coords[i]);
//          Vector2D pt_vec = Vector2D(x, y);
//          crosses[i] = edge_vec.x * pt_vec.y - pt_vec.y * edge_vec.x;
//
//
//      }
//      std::cout << crosses << std::endl;
//
//      if ((crosses[0] < 0 && crosses[1] < 0 && crosses[2] < 0) || (crosses[0] > 0 && crosses[1] > 0 && crosses[2] > 0)) {
//          return 1;
//      }
//      return 0;
//A: 0, B: 1
//    Vector2D edges[] = {Vector2D(x0-x1,y0-y1), Vector2D(x0-x2,y0-y2), Vector2D(x1-x2, y1-y2)};
    Vector2D l_norm[] = {Vector2D(y0-y1, x1-x0), Vector2D(y0-y2, x2-x0), Vector2D(y1-y2, x2-x1)};
    bool line_check[3];

    for (int i = 0; i < 3; i++) {
        Vector2D norm = l_norm[i];
        Vector2D pt = Vector2D(abs(x_coords[i]-x), y_coords[i]-y);
        line_check[i] = norm.x * pt.x + norm.y + pt.y;
        if ((norm.x * pt.x + norm.y + pt.y) >= 0) {
            return 0;
        }
    }
    return 1;

  }

  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    int xmin = floor(min(min(x0, x1), x2));
    int xmax = floor(max(max(x0, x1), x2));
    int ymin = floor(min(min(y0, y1), y2));
    int ymax = floor(max(max(y0, y1), y2));


    for (int x = xmin; x <= xmax; x++) {
        for (int y = ymin; y <= ymax; y++) {
            if (is_in_triangle(x0, y0,x1, y1, x2, y2,x+0.5,y+0.5)) {
                rasterize_point(x+0.5, y+0.5, color);
            }
        }
    }


    // TODO: Task 2: Update to implement super-sampled rasterization


  return;


  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle



  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
