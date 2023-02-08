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

    sample_buffer[y * width * sqrt(sample_rate) + x] = c;

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

    // fill in supersampled buffer
    float sqrt_rate = sqrt(sample_rate);
    for (float i = 0; i < sqrt_rate; i++) {
        for (float j = 0; j < sqrt_rate; j++) {
            fill_pixel((sx * sqrt_rate + i), (sy * sqrt_rate + j), color);
        }
    }

    //fill_pixel(sx, sy, color);
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

  // Created for Task 1
  float L(float x, float y, float P1x, float P1y, float P2x, float P2y) {
      return (x - P2x) * (P1y - P2y) - (y - P2y) * (P1x - P2x);
  }

  // Created for Task 1
  bool is_inside_triangle(float P1x, float P1y, float P2x, float P2y, float P3x, float P3y, float x, float y) {
      float l1, l2, l3;
      l1 = L(x, y, P1x, P1y, P2x, P2y);
      l2 = L(x, y, P2x, P2y, P3x, P3y);
      l3 = L(x, y, P3x, P3y, P1x, P1y);
      return ((l1 >= 0) && (l2 >= 0) && (l3 >= 0)) || ((l1 < 0) && (l2 < 0) && (l3 < 0));
  }

  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {

    // Min and max edges for trying sample points
    int xmin = floor(min(min(x0, x1), x2));
    int xmax = floor(max(max(x0, x1), x2));
    int ymin = floor(min(min(y0, y1), y2));
    int ymax = floor(max(max(y0, y1), y2));

    // Sample rate
    float sqrt_rate = sqrt(sample_rate);
    float sample_size = 1.0 / sqrt_rate;
    float center = sample_size / 2.0;

    // Iterate through pixels and rasterize
    for (int x = xmin; x <= xmax; x++) {
        for (int y = ymin; y <= ymax; y++) {
            for (float i = 0; i < sqrt_rate; i++) {
                for (float j = 0; j < sqrt_rate; j++) {
                    float px = x + (i * sample_size) + center;
                    float py = y + (j * sample_size) + center;
                    if (is_inside_triangle(x0, y0, x1, y1, x2, y2, px, py)) {
                        fill_pixel((x * sqrt_rate + i), (y * sqrt_rate + j), color);
                    }
                }
            }
        }
    }


    // Done: Task 2: Update to implement super-sampled rasterization


  return;


  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // Done: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    // Min and max edges for trying sample points
      int xmin = floor(min(min(x0, x1), x2));
      int xmax = floor(max(max(x0, x1), x2));
      int ymin = floor(min(min(y0, y1), y2));
      int ymax = floor(max(max(y0, y1), y2));

      // Sample rate
      float sqrt_rate = sqrt(sample_rate);
      float sample_size = 1.0 / sqrt_rate;
      float center = sample_size / 2.0;

      // Iterate through pixels and rasterize
      for (int x = xmin; x <= xmax; x++) {
          for (int y = ymin; y <= ymax; y++) {
              for (float i = 0; i < sqrt_rate; i++) {
                  for (float j = 0; j < sqrt_rate; j++) {
                      float px = x + (i * sample_size) + center;
                      float py = y + (j * sample_size) + center;
                      if (is_inside_triangle(x0, y0, x1, y1, x2, y2, px, py)) {
                          float alpha, beta, gamma;
                          alpha = ((py - y1) * (x2 - x1) - (px - x1) * (y2 - y1)) / ((y0 - y1) * (x2 - x1) - (x0 - x1) * (y2 - y1));
                          beta = ((py - y2) * (x0 - x2) - (px - x2) * (y0 - y2)) / ((y1 - y2) * (x0 - x2) - (x1 - x2) * (y0 - y2));
                          gamma = 1 - alpha - beta;
                          Color c = alpha * c0 + beta * c1 + gamma * c2;
                          fill_pixel((x * sqrt_rate + i), (y * sqrt_rate + j), c);
                      }
                  }
              }
          }
      }

      return;

  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // Done: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    // Min and max edges for trying sample points
      int xmin = floor(min(min(x0, x1), x2));
      int xmax = floor(max(max(x0, x1), x2));
      int ymin = floor(min(min(y0, y1), y2));
      int ymax = floor(max(max(y0, y1), y2));

      // Sample rate
      float sqrt_rate = sqrt(sample_rate);
      float sample_size = 1.0 / sqrt_rate;
      float center = sample_size / 2.0;

      // Iterate through pixels and rasterize
      for (int x = xmin; x <= xmax; x++) {
          for (int y = ymin; y <= ymax; y++) {
              for (float i = 0; i < sqrt_rate; i++) {
                  for (float j = 0; j < sqrt_rate; j++) {
                      float px = x + (i * sample_size) + center;
                      float py = y + (j * sample_size) + center;
                      if (is_inside_triangle(x0, y0, x1, y1, x2, y2, px, py)) {
                          float alpha, beta, gamma;
                          alpha = ((py - y1) * (x2 - x1) - (px - x1) * (y2 - y1)) / ((y0 - y1) * (x2 - x1) - (x0 - x1) * (y2 - y1));
                          beta = ((py - y2) * (x0 - x2) - (px - x2) * (y0 - y2)) / ((y1 - y2) * (x0 - x2) - (x1 - x2) * (y0 - y2));
                          gamma = 1 - alpha - beta;

                          SampleParams sp;
                          sp.p_uv = Vector2D(alpha * u0 + beta * u1 + gamma * u2, alpha * v0 + beta * v1 + gamma * v2);
                          sp.p_dx_uv = Vector2D((x0 - x1)/(u0 - u1), (x0 - x1) / (v0 - v1));
                          sp.p_dy_uv = Vector2D((y0 - y1) / (u0 - u1), (y0 - y1) / (v0 - v1));
                          sp.psm = psm;
                          sp.lsm = lsm;

                          Color c;
                          c = tex.sample(sp);
                          fill_pixel((x * sqrt_rate + i), (y * sqrt_rate + j), c);
                      }
                  }
              }
          }
      }


  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // Done: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * rate, Color::White);
    clear_buffers();
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // Done: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    //this->sample_buffer.resize(width * height, Color::White);
    set_sample_rate(this->sample_rate);
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
    // Done: Task 2: You will likely want to update this function for supersampling support


     float rate = sqrt(sample_rate);

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {

          // Average high-res sample values
          Color col = Color::Black;
          for (int i = 0; i < rate; i++) {
              for (int j = 0; j < rate; j++) {
                  col += sample_buffer[(rate * y + j) * width * sqrt(sample_rate) + (rate * x + i)];
              }
          }
          col *= (1.0 / sample_rate);
        

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
