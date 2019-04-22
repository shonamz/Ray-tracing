/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Namazian>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <string>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define PI 3.14159265

const double BIS= 1e-16;

const double MDIST = -1e8;

const int MAXIMUM_SCREENSHOTS = 10;

int frameNum;

int g_screenshotCounter = 0;

//screen shot

bool takeScreenshots = true;


unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};



Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Color definition
struct Color {
  double vR;
  double vG;
  double vB;

  Color () : vR(0), vG(0), vB(0) {}
  Color (double r, double g, double b) : vR(r), vG(g), vB(b) {}

  Color& operator += (const Color& other) {
    // Clamp color  between 0 and 1
    vR += other.vR;
    if (vR > 1.0f) {
      vR = 1.0f;
    }

    else if (vR < 0) {
      vR = 0;
    }

    vG += other.vG;
    if (vG > 1.0f) {
      vG = 1.0f;
    }

    else if (vG < 0) {
      vG = 0;
    }

    vB += other.vB;
    if (vB > 1.0f) {
      vB = 1.0f;
    }

    else if (vB < 0) {
      vB = 0;
    }

    return *this;
  }
};
 


struct Vector3{

  double vX;
  double vY;
  double vZ;

  Vector3 () : vX(0), vY(0), vZ(0) {}
  Vector3 (double x, double y, double z) : vX(x), vY(y), vZ(z) {}

    // function 
  static Vector3  crossProduct (const Vector3& a, const Vector3& b);

  inline double distance () const { return std::pow(vX, 2) + std::pow(vY, 2) + std::pow(vZ, 2); }
  inline double dotProduct(const Vector3 & un ) const { return  vX * un.vX + vY * un.vY +vZ * un.vZ;}
  Vector3& normalize();
inline double magnitude () const { return std::sqrt(distance()); }


  //operator
  Vector3 operator- (const Vector3& un) const { return Vector3 (vX - un.vX, vY - un.vY, vZ - un.vZ); }
  Vector3 operator* (double scal) const { return Vector3 (vX * scal, vY * scal, vZ * scal); }
  Vector3 operator+ (const Vector3& un) const { return Vector3 (vX + un.vX, vY + un.vY, vZ + un.vZ); }
  Vector3 operator- () const { return Vector3 (-vX, -vY, -vZ); }



  Vector3& operator += (const Vector3& un) {
    vX += un.vX;
    vY += un.vY;
    vZ += un.vZ;
    return *this;
  }

 };

Vector3& Vector3::normalize () {
 //check to not divided by zero
  double norm = distance();
  if (norm > 0) {
    double inver = 1.0f / std::sqrt(norm);
    vX *= inver;
    vY *= inver;
    vZ *= inver;
  }
  return *this;
}


Vector3 Vector3::crossProduct(const Vector3& a ,const Vector3& b ){

  double x = a.vY * b.vZ - a.vZ * b.vY;
  double y = a.vZ * b.vX - a.vX * b.vZ;
  double z = a.vX * b.vY - a.vY * b.vX;

  return Vector3 (x, y, z);

}

class Ray {

 private :

 Vector3 vDirection;
 Vector3 vOrigin;

  public:

  Ray () {}
  Ray (const Vector3& origin, const Vector3& direction) : vOrigin (origin), vDirection (direction) {}

  // Member functions
  bool intersect (const Sphere& sphere, Vector3& intersection);
  bool intersect (const Triangle& triangle, Vector3& intersection);

};

bool Ray::intersect(const Sphere& sphere, Vector3& intersection){

// center of the sphere and the distance vector  

  Vector3 position = Vector3(sphere.position[0], sphere.position[1], sphere.position[2]);
  Vector3 dist = vOrigin - position;

  // Create values that correspond to the results of the quadratic equation. There are two.

  double t0= MDIST;
  double t1=MDIST ;

  // Calculate the quadratic equations for the top and bottom of the sphere
  double a = vDirection.dotProduct(vDirection);
  double b = 2 * vDirection.dotProduct (dist);
  double c = dist.dotProduct (dist) - std::pow (sphere.radius, 2);
  double quad = std::pow(b, 2) - (4 * a * c);
  if (quad < 0) {
    return false;
  }

  //  discriminant is nearly close  to 0, t0 = t1 = -0.5 * b/a;

  else if (std::abs(quad) < BIS) {
    double q = -0.5f * b / a;
    t0 = q;
    t1 = q;
  }

  //  discriminant is greater than 0
  else {
    double r = (b > 0) ? -0.5f * (b + std::sqrt(quad)) : -0.5f * (b - std::sqrt(quad));
    t0 = r / a;
    t1 = c / r; 
  }

  // Check  t values are positive or 0

  if (t0 < 0 && t1 < 0) {
    return false;
  }
// return min (t0,t1)

  if (t0 > t1 && t1 > 0) {
    t0 = t1;
  }

  intersection = vOrigin + (vDirection * t0);
  return true;
}

bool Ray::intersect (const Triangle& triangle, Vector3& intersection){

  // 
 Vector3 verA= Vector3(triangle.v[0].position[0],triangle.v[0].position[1],triangle.v[0].position[2]);
 Vector3 verB= Vector3(triangle.v[1].position[0],triangle.v[1].position[1],triangle.v[1].position[2]);
 Vector3 verC= Vector3(triangle.v[2].position[0],triangle.v[2].position[1],triangle.v[2].position[2]);


//normal for ray-plane intersection 

Vector3 normalPlane = Vector3::crossProduct(verB - verA, verC - verA);
normalPlane.normalize ();


// do not product to chekck if plane and ray are parallel  dot product is zero or close to zero 

  double direction = normalPlane.dotProduct(vDirection);
  if (std::abs(direction) < BIS) {
    return false;
  }

// d parameter

  Vector3 distFromOrigin = verA - vOrigin;
  double d = distFromOrigin.dotProduct(normalPlane);

//  h value is negative (behind the ray) no intersection

  double t = d / direction;
  if (t < 0) {
    return false;
  }
  // intersection point

  intersection = vOrigin + (vDirection * t);

  // This is an in out test

  if (normalPlane.dotProduct(Vector3::crossProduct(verB - verA, intersection - verA)) < 0
    || normalPlane.dotProduct(Vector3::crossProduct(verC - verB, intersection - verB)) < 0
    || normalPlane.dotProduct(Vector3::crossProduct(verA - verC, intersection - verC)) < 0) {
    return false;
  }

  return true;

}

// Find ight direction dot normal for traingle and sphere
double computeLight (const Vector3& lightDir, const Vector3& normal) {

  // Find the dot product of the light direction and the normal L dot N  then clamp it
  double lightdNorm = lightDir.dotProduct(normal);
  if ( lightdNorm> 1.0f) {
    lightdNorm = 1.0f;
  }

  else if (lightdNorm < 0.0f) {
    lightdNorm = 0.0f;
  }

  return lightdNorm;
}
// Find the magnitude for  reflectiveness for spheres and triangles
double computeReflect (double lightMag, const Vector3& lightDir, const Vector3& intersection, const Vector3& normal) {

  // Get normalized direction vector
  Vector3 direction = -intersection;
  direction.normalize ();

  // Find the reflection vector, then compute the dot product to get the reflection magnitude
  Vector3 reflection (2 * lightMag * normal.vX - lightDir.vX, 2 * lightMag * normal.vY - lightDir.vY, 2 * lightMag * normal.vZ - lightDir.vZ);
  reflection.normalize ();
  double rDotV = reflection.dotProduct(direction);
  if (rDotV > 1.0f) { 
    rDotV = 1.0f;
  }

  else if (rDotV < 0.0f) {
    rDotV = 0.0f;
  }

  return rDotV;
}


Color sphereLighting (const Sphere& sphere, const Light& light, const Vector3& intersection) {

  // normal
  Vector3 normal = intersection - Vector3 (sphere.position[0], sphere.position[1], sphere.position[2]);
  normal.normalize ();

  //  normalized light direction 
  Vector3 lightPosition (light.position[0], light.position[1], light.position[2]);
  Vector3 lightDir = lightPosition - intersection;
  lightDir.normalize ();

  // Compute and clamp LdotN  and (2 * l (n - dir)) * direction (the reflection magnitude)
  double lightMagnitude = computeLight(lightDir, normal);
  double reflection = computeReflect(lightMagnitude, lightDir, intersection, normal);

  // Get diffuse, specular, and shininess  from  sphere object
  Color diffuse (sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
  Color specular (sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
  double shinin = sphere.shininess;

  //  intensity for each color by phong 
  double r = light.color[0] * (diffuse.vR * lightMagnitude + (specular.vR * std::pow (reflection, shinin)));
  double g = light.color[1] * (diffuse.vG * lightMagnitude + (specular.vG * std::pow (reflection, shinin)));
  double b = light.color[2] * (diffuse.vB * lightMagnitude + (specular.vB * std::pow (reflection, shinin)));

  return Color (r, g, b);
}


Color triangleLighting (const Triangle& triangle, const Light& light, const Vector3& intersection) {

  // To find normal, the barycentric coordinates need to be computed
  Vector3 vertexA (triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
  Vector3 vertexB (triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
  Vector3 vertexC (triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);


  // Find planar normal
  Vector3 aB = vertexB - vertexA;
  Vector3 aC = vertexC - vertexA;

  Vector3 pla = Vector3::crossProduct (aB, aC);
  float denom = pla.dotProduct(pla);

  // Find values u and v
  Vector3 edgeCB = vertexC - vertexB;
  Vector3 distB = intersection - vertexB;
  Vector3 cpCB = Vector3::crossProduct(edgeCB, distB);

  Vector3 edgeAC = vertexA - vertexC;
  Vector3 distC = intersection - vertexC;
  Vector3 cpAC = Vector3::crossProduct(edgeAC, distC);

  // Find final u, v, and w values
  double u = pla.dotProduct(cpCB) / denom;  // Alpha
  double v = pla.dotProduct(cpAC) / denom;  // Beta
  double w = 1.0f - u - v;          // Gamma

  // Find triangle normals
  Vector3 normal (
    u * triangle.v[0].normal[0] + v * triangle.v[1].normal[0] + w * triangle.v[2].normal[0],
    u * triangle.v[0].normal[1] + v * triangle.v[1].normal[1] + w * triangle.v[2].normal[1],
    u * triangle.v[0].normal[2] + v * triangle.v[1].normal[2] + w * triangle.v[2].normal[2]
  );
  normal.normalize ();

  // Find normalized light direction vector
  Vector3 lightPosition (light.position[0], light.position[1], light.position[2]);
  Vector3 lightDirection = lightPosition - intersection;
  lightDirection.normalize ();

  // Find  and clamp the values of the magnitudes of LdotN (light magnitude) and (2 * l (n - dir)) * direction (the reflection magnitude)
  double lightMagnitude = computeLight (lightDirection, normal);
  double reflectionMagnitude = computeReflect(lightMagnitude, lightDirection, intersection, normal);

  // Find diffuse, specular, and shininess values from the triangle object
  Color diffuse (
    u * triangle.v[0].color_diffuse[0] + v * triangle.v[1].color_diffuse[0] + w * triangle.v[2].color_diffuse[0],
    u * triangle.v[0].color_diffuse[1] + v * triangle.v[1].color_diffuse[1] + w * triangle.v[2].color_diffuse[1],
    u * triangle.v[0].color_diffuse[2] + v * triangle.v[1].color_diffuse[2] + w * triangle.v[2].color_diffuse[2]
  );

  Color specular (
    u * triangle.v[0].color_specular[0] + v * triangle.v[1].color_specular[0] + w * triangle.v[2].color_specular[0],
    u * triangle.v[0].color_specular[1] + v * triangle.v[1].color_specular[1] + w * triangle.v[2].color_specular[1],
    u * triangle.v[0].color_specular[2] + v * triangle.v[1].color_specular[2] + w * triangle.v[2].color_specular[2]
  );

  double shininess = u * triangle.v[0].shininess + v * triangle.v[1].shininess + w * triangle.v[2].shininess;

  // Find intensity for each color using the Phong equation
  double r = light.color[0] * (diffuse.vR * lightMagnitude + (specular.vR * std::pow (reflectionMagnitude, shininess)));
  double g = light.color[1] * (diffuse.vG * lightMagnitude + (specular.vG * std::pow (reflectionMagnitude, shininess)));
  double b = light.color[2] * (diffuse.vB * lightMagnitude + (specular.vB * std::pow (reflectionMagnitude, shininess)));
  
  return Color (r, g, b);
}

// Go through all spheres and set the pixel color depending base on  intersection and shadowing
Color sphereTest (Ray& ray, const Color& color, double& closest) {
 
  Color ret = color;
  for (int i = 0; i < num_spheres; i++) {

    // find  color for the closest intersection and overwrite ret if a closer intersection is detected
    Vector3 intersection (0, 0, MDIST);
    if (ray.intersect (spheres[i], intersection) && intersection.vZ > closest) {

      // By default, the color should be black
      ret = Color (0, 0, 0);

      // Check to see if the sphere is shadowed--if it isn't, add the color at the intersection point
      for (int j = 0; j < num_lights; j++) {

        // Get the position of the light
        Vector3 lightPosition (lights[j].position[0], lights[j].position[1], lights[j].position[2]);

        // Create the shadow ray--the origin should be the point where the ray intersected with the sphere, and the direction should be the normalized direction to the light
        Vector3 origin = intersection;
        Vector3 direction = lightPosition - origin;
        Ray shadow (origin, direction.normalize ());

        // If the object is lit we are going to add color to it
        bool isOn = true;

        // Check for collisions against every object
        for (int k = 0; k < num_spheres; k++) {
          Vector3 intersect;
          if (shadow.intersect (spheres[k], intersect) && k != i) {

            // check that the intersection  is not past the light
            Vector3 x = intersect - intersection;
            Vector3 y = lightPosition - intersection;
            if (x.magnitude () < y.magnitude ()) {
              isOn = false;
              break;
            }
          }
        }

        for (int k = 0; k < num_triangles; k++) {
          Vector3 intersect;
          if (shadow.intersect(triangles[k], intersect)) {

            // check if the intersection point not past the light
            Vector3 x = intersect - intersection;
            Vector3 y = lightPosition - intersection;
            if (x.magnitude () < y.magnitude ()) {
              isOn = false;
              break;
            }
          }
        }

        if (isOn) {
          ret += sphereLighting (spheres[i], lights[j], intersection);
 

        }
      }

      // Update the closest intersection
      closest = intersection.vZ;
    }
  }
 
 
  return ret;
}


// go  through all triangles to set the pixel color
Color triangleTest (Ray& ray, const Color& color, double& closest) {
 
  Color ret = color;
  for (int i = 0; i < num_triangles; i++) {

    //   color for the closest intersection but overwrite ret if  find a closer intersection
    Vector3 intersection (0, 0, MDIST);
    if (ray.intersect(triangles[i], intersection) && intersection.vZ > closest) {

      // default the color is black
      ret = Color (0, 0, 0);

      // if the triangle is in  shadow-- isn'if not  add the color at the intersection
      for (int j = 0; j < num_lights; j++) {

        // Get the position of the light
        Vector3 lightPosition (lights[j].position[0], lights[j].position[1], lights[j].position[2]);

        // Create the shadow ray--the origin should be the point where the ray intersected with the sphere, and the direction should be the normalized direction to the light
        Vector3 origin = intersection;
        Vector3 direction = lightPosition - origin;
        Ray shadow (origin, direction.normalize ());

        //  add color to it if it is lit
        bool isOn = true;

        // Check for collisions against every object--ignoring our own object (same index)
        for (int k = 0; k < num_spheres; k++) {
          Vector3 intersect;
          if (shadow.intersect(spheres[k], intersect)) {

            // Make sure that the intersection point is not past the light
            Vector3 x = intersect - intersection;
            Vector3 y = lightPosition - intersection;
            if (x.magnitude () < y.magnitude ()) {
              isOn = false;
              break;
            }
          }
        }

        for (int k = 0; k < num_triangles; k++) {
          Vector3 intersect;
          if (shadow.intersect(triangles[k], intersect) && k != i) {

            // Make sure that the intersection point is not past the light
            Vector3 x = intersect - intersection;
            Vector3 y = lightPosition - intersection;
            if (x.magnitude () < y.magnitude ()) {
              isOn = false;
              break;
            }
          }
        }

        if (isOn) {
          ret += triangleLighting (triangles[i], lights[j], intersection);
        }
      }

      // Update the closest intersection
      closest = intersection.vZ;

    }
  }
 

  return ret;
}


//  raycast to find a pixel's color
Color cast (Ray& ray) {
  
  // find the current pixel color and current closest intersection
  Color retCast (1, 1, 1);

  double closest = MDIST;

  // raycasts for all spheres
  retCast = sphereTest (ray, retCast, closest);

  // raycasts for all triangles
  retCast = triangleTest (ray, retCast, closest);

  // Add ambient light
  retCast += Color (ambient_light[0], ambient_light[1], ambient_light[2]);

 
  // in case  no triangles or spheres, return
  return retCast;
}



Ray generateRayFromCamera(double x, double y){

  // the ratio and angle

  double ratio = (double)WIDTH / (double)HEIGHT;
  double angle = std::tan ((PI / 180.0) * (fov / 2.0));

// NDC coordinates

  double xN = (x + 0.5) / (double)WIDTH;
  double yN = (y + 0.5) / (double)HEIGHT;

  double xS = (2 * xN) - 1;
  double yS = (2 * yN) - 1;

  // with the screen coordinates, get the real  coordinates

  double xReal = xS * angle * ratio;
  double yReal = yS * angle;

  Vector3 origin;
  Vector3 direction (xReal, yReal, -1);

   return  Ray (origin, direction.normalize () );

    
}

// FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {

// calculate four rays
      double r = 0;
      double g = 0;
      double b = 0;

        Ray ray = generateRayFromCamera(x, y);
        Color color = cast (ray);
        r = color.vR;
        g = color.vG;
        b = color.vB;
       plot_pixel(x, y, r * 255, g * 255, b * 255);

    }
    glEnd();
    glFlush();
  }
  printf("Done!\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

// write a screenshot to the specified filename
void save_jpg(const char * filename)
{
    unsigned char * screenshotData = new unsigned char[WIDTH * HEIGHT * 3];
    
    glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, screenshotData);
    
    ImageIO img(WIDTH, HEIGHT, 3, screenshotData);

   // ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
    
    
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}



void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
    //  save_jpg();
      if (g_screenshotCounter < MAXIMUM_SCREENSHOTS && takeScreenshots) {
          std::string filename = std::to_string(frameNum) + ".jpg";
          if (frameNum < 100) {
              filename = "0" + filename;
          }
          if (frameNum < 10) {
              filename = "0" + filename;
          }
          filename = "Recording/" + filename;
          save_jpg(filename.c_str());
          ++frameNum;
          g_screenshotCounter++;
      }
  }
// once=1;
    
    


}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

