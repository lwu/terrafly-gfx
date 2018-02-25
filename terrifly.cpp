//
// Terrifly!
// 
// Leslie Wu
//
//

#include <gfx/gui.h>
#include <gfx/raster.h>
#include <gfx/mat4.h>
#include <gfx/geom3d.h>

#include <FL/fl_file_chooser.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/filename.H>
#include <FL/fl_ask.H>

#include <mmsystem.h>
#include <assert.h>
#include <deque>

#include "precalc.h"
#include "terrain.h"
#include "error.h"

#include "globals.h"
#include "glwin.h"

// global options
bool enable_lighting = false;
bool display_error_lines = false;
const float VIEW_ANGLE = 60.0f;
bool recalc_active = true;
bool timedemo = false;
bool use_vertex_array = false;
bool ground_walk = false;
bool ground_collide = true;
bool fly = false;
bool adaptive_tau = false;
bool guarantee_tri_limit = false;
bool gps = false;
bool rot_compensate = true;
bool phi_compensate = false;

using namespace std;


extern void error_popup(const char* err_msg)
{
	fl_alert(err_msg);
	exit(0);
}


Vec3f rotate_direction(const Vec3f& moving, float angle, const Vec3f& axis)
{
	Vec3 daxis(axis[0], axis[1], axis[2]);
	Vec3 dmoving(moving[0], moving[1], moving[2]);
    Mat4 M = rotation_matrix_deg(angle, daxis);
    Vec4 v = Vec4(dmoving, 0.0);
    Vec4 vnew = M*v;

    return Vec3f(vnew[0], vnew[1], vnew[2]);
}




class GUI : public MxGUI
{

public:
    Terrain terrain;

    bool will_draw_wireframe;
    bool will_draw_texture;

	bool update_terrain_viewpoint;

    float move_step;

    float znear, zfar;
    Vec3f eye, forward, up;
	Vec3f forward_2d;
	float theta, phi;
	float dtheta;
	bool mouse_move_forward, mouse_pressed;
	int mouse_x0, mouse_y0;
	int mouse_x1, mouse_y1;

    GUI();

    virtual void setup_for_drawing();
    virtual void draw_contents();
    virtual bool key_press(int key);
    virtual void cmdline_file(const char *filename);

	virtual bool mouse_drag(int *where, int *last, int which);
	virtual bool mouse_down(int *where, int which);
	virtual bool mouse_up(int *where, int which);
	virtual void update_animation();

	virtual void add_upper_controls(int& yfill, const int pad);

    void apply_camera();

	void update_rotation();
	void update_pos();

public:
	Fl_Value_Slider* tau_slider;
	Fl_Value_Slider* tri_slider;
	Fl_Value_Slider* height_slider;
	Fl_Value_Slider* speed_slider;
};

GUI gui;

GUI::GUI()
{
    znear = 5.0f;
    zfar = 6000.0f;
    will_draw_wireframe = true;
    will_draw_texture = false;

	theta = -63.44f;
	phi = 0.0f;
	dtheta = 0.0f;
	mouse_move_forward = false;

	update_terrain_viewpoint = true;
}



void GUI::setup_for_drawing()
{	
	float rgb[] = {37.0f, 34.0f, 172.0f};
    glClearColor(rgb[0]/255.0f, rgb[1]/255.0f, rgb[2]/255.0f, 0.0f);

    glEnable(GL_DEPTH_TEST);

	if (enable_lighting) {
		glEnable(GL_NORMALIZE);
		
		// Enable lighting and set up the lighting environment.  We specify a
		// global ambient glow and create two point lights.
		
		glEnable(GL_LIGHTING);
		Vec4f ambient_light(1.0, 1.0, 1.0, 1.0);
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (float *)ambient_light);
		
		const Vec4f light0_pos(0.0f, 0.5f, 1.0f, 0.0f);
		glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
		glEnable(GL_LIGHT0);
		
		const Vec4f rgb(0.912f, 0.717f, 0.505f, 1.0f);
		const Vec4f r_amb  = 0.1*rgb;
		const Vec4f r_diff = 1.0*rgb;
		const Vec4f r_spec = 0.3*rgb;
		
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, r_amb);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, r_diff);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, r_spec);
		glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 100);
	}

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	apply_camera(); // setup some vars

	static bool inited = false;

	if (inited) return;

	setup_gl_extensions();

    if( terrain.texture )
    {
	// Can either use the texture to MODULATE the lighting of the
	// underlying surface, or just entirely REPLACE lighting
	// calculations with the color in the texture map.
	//
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

#if !defined(NO_MIPMAP)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);	
	
	GLint internalFmt = GL_RGB;
	gluBuild2DMipmaps(GL_TEXTURE_2D,
			  internalFmt,
			  terrain.texture->width(),
			  terrain.texture->height(),
			  GL_RGB,			  
			  GL_UNSIGNED_BYTE,
			  terrain.texture->head());
#else	
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	
	
	// Without mip-mapping use:
	
	  glTexImage2D(GL_TEXTURE_2D, 0, terrain.texture->channels(),
	               terrain.texture->width(), terrain.texture->height(),
		       0, GL_RGB, GL_UNSIGNED_BYTE,
		       terrain.texture->head());
#endif	

	delete terrain.texture;	

    }

	inited = true;
}

void GUI::apply_camera()
{
    float aspect = (float)canvas->w() / (float)canvas->h();

    glMatrixMode(GL_PROJECTION);
    gluPerspective(VIEW_ANGLE, aspect, znear, zfar);

    glMatrixMode(GL_MODELVIEW);

    Vec3f at = eye + forward;
    gluLookAt(eye[0], eye[1], eye[2],
              at[0], at[1], at[2],
              up[0], up[1], up[2]);

	terrain.draw_style_fast = !guarantee_tri_limit;

	if (update_terrain_viewpoint) {
		terrain.set_eye(eye, forward_2d, up, VIEW_ANGLE, aspect, dtheta, phi);
		terrain.tri_limit = (int)(tri_slider->value() * 1000.0f);
	}	

	float new_tau = tau_slider->value();

	// Simple control system that increases or decreases tau (pixel error)
	//  to get # of triangles drawn closer to triangle limit
	if (adaptive_tau) {
		int delta_tri = terrain.total_tri_count - terrain.tri_limit;

		if (delta_tri >= -500) {
			// surplus triangles
			new_tau += 0.02f;
			new_tau = min(16.0f, new_tau);
		} else if (delta_tri < -(terrain.tri_limit/2)) {
			// not enough triangles
			new_tau -= 0.01f;
			new_tau = max(1.0f, new_tau);
		}
	}

	tau_slider->value(new_tau);
	int new_lambda = canvas->w();
	terrain.set_error_tolerance(new_lambda, new_tau);
}


// GPS system in upper-right hand corner
void draw_gps()
{		
	glDisable(GL_DEPTH_TEST);
	float w = (float)gui.canvas->w();
	float h = (float)gui.canvas->h();

	glMatrixMode(GL_PROJECTION); glLoadIdentity();
	gluOrtho2D(0.0, 1.0f, 0.0f, 1.0f);
	glMatrixMode(GL_MODELVIEW); glLoadIdentity();
	
	float gps_width  = 0.15f;
	float gps_height = 0.12f;

	float gps_x = 0.80f;
	float gps_y = 0.80f;

	glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);
		glVertex2f(gps_x, gps_y);
		glTexCoord2f(1.0f, 0.0f);
		glVertex2f(gps_x + gps_width, gps_y);
		glTexCoord2f(1.0f, 1.0f);
		glVertex2f(gps_x + gps_width, gps_y + gps_height);
		glTexCoord2f(0.0f, 1.0f);
		glVertex2f(gps_x, gps_y + gps_height);
	glEnd();

	float x = gui.eye[0] * gui.terrain.inv_w / gui.terrain.pixel_spacing;
	float y = gui.eye[1] * gui.terrain.inv_h / gui.terrain.pixel_spacing;

	if ( gui.will_draw_texture ) glDisable(GL_TEXTURE_2D);    

	glColor3f(0.0f, 1.0f, 0.0f);
	glPointSize(3.0f);
	glBegin(GL_POINTS);
		glVertex2f(gps_x + x*gps_width, gps_y + y*gps_height);
	glEnd();
	glColor3f(1.0f, 1.0f, 1.0f);

	if ( gui.will_draw_texture ) glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);
}

// class used to average framerate. whee.
template <typename T, int N>
struct n_queue {	
	n_queue() {
		cur_index = 0;
		for (int i=0; i < N; i++) {
			que[i] = T();
		}
	}

	void add(const T& elem) {
		que[cur_index] = elem;
		cur_index = (cur_index+1) % N;
	}

	T avg() const {		
		T val = T();
		for (int i=0; i < N; i++) {
			val += que[i];			
		}

		return val / (T)N;
	}

private:	
	int cur_index;
	T que[N];	
};

n_queue<float, 32> fps_queue;

void GUI::draw_contents()
{
	int loop_count = timedemo ? 90 : 1;
	
	double t_0 = get_cpu_time();

	// main terrain rendering loop. to measure performance, 
	//  we loop multiple times and record the average fps
	for (int i=0; i<loop_count; i++) {		
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
		
		if( will_draw_texture )						glEnable(GL_TEXTURE_2D);
		else                                        glDisable(GL_TEXTURE_2D);
		
		if( will_draw_wireframe ) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		else                      glPolygonMode(GL_FRONT, GL_FILL);
		//else                      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				
		glMatrixMode(GL_PROJECTION);  glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);   glLoadIdentity();
		
		apply_camera();
		
		terrain.draw();
		
		if (timedemo) {
			eye += forward * 32.0f;
		}
		
		if (gps) {		
			draw_gps();
		}
	}

	double t_1 = get_cpu_time();		
	double dt = t_1 - t_0;
	float fps;	

	if (dt > 0) {		
		fps = (float)loop_count / (float)dt;
	} else {
		fps = 120.0f;
	}

	fps_queue.add(fps);
	float avg_fps = fps_queue.avg();

	if (timedemo) {
		status("timedemo: %2.1f average fps", fps);		
		timedemo = false;
	} else {
		status("%05.1f fps, %05.1f avg fps, subdiv lvl %d, %d tris, %d updated tris", 
			fps, avg_fps, terrain.subdiv_level, terrain.total_tri_count, terrain.tri_count
		);
	}	
}

void normalize(Vec3f& v)
{
	v /= norm(v);
}

bool GUI::mouse_down(int *where, int which)
{
	mouse_x0 = where[0];
	mouse_y0 = where[1];
	
	mouse_pressed = true;

	if (which > 1) {
		mouse_move_forward = true;
	}
	return true;
}

bool GUI::mouse_up(int *where, int which)
{
	if (which > 1) {
		mouse_move_forward = false;
	}

	mouse_pressed = false;

	return true;
}

void GUI::update_rotation() {
	Vec3f dir(0.0f, 1.0f, 0.0f);
	forward = rotate_direction(dir, theta, up);

	forward_2d = forward; // keep 2-D forward for view culling

	Vec3f left = cross(up, forward);
	forward = rotate_direction(forward, phi, left);	
}

void GUI::update_pos()
{	
	// drop camera to ground if option is enabled
	if (ground_collide || ground_walk) {
		int x = eye[0] * 1.0f / terrain.pixel_spacing;
		int y = eye[1] * 1.0f / terrain.pixel_spacing;

		static float dy = 0.0f;

		if (x >= 0 && x < terrain.z->width() &&
			y >= 0 && y < terrain.z->height()) {
			float ground_z = terrain.z->pixel(x,y)[0] * terrain.height_factor;			

			float new_eye_z = ground_z + height_slider->value();

			if (ground_walk) { //gravity
				if (eye[2] > new_eye_z) {
					dy += -1.0f;
					eye[2] += dy;
				}
			} 

			if (eye[2] < new_eye_z) {
				eye[2] = new_eye_z;			
				dy = 0.0f;
			}			
		}
	}

	// maintain min distance
	if (fly) {
		eye[2] = max(eye[2], height_slider->value());
	}

}

bool GUI::mouse_drag(int *where, int *last, int which)
{
	mouse_x1 = where[0]; mouse_y1 = where[1];

	return true;
}

void bound_range(float* x, float x_low, float x_high)
{
	if (*x < x_low) {
		*x = x_low;
	} else if (*x > x_high) {
		*x = x_high;
	}
}

void GUI::update_animation()
{
	// update camera motion if mouse is active
	if (mouse_move_forward) {
		move_step = speed_slider->value();

		eye += move_step * forward;
		update_pos();
	}

	if (mouse_pressed) {
		int x0 = mouse_x0, y0 = mouse_y0;
		int x1 = mouse_x1, y1 = mouse_y1;
		int dx = x1 - x0, dy = y1 - y0;
		
		dtheta = -(dx * 0.01f);

		theta += dtheta;
		phi += dy * 0.01f;

		bound_range(&phi, -85.0f, 60.0f);
		
		update_rotation();
	} else {
		dtheta = 0.0f;
	}
}

bool GUI::key_press(int key)
{
	move_step = speed_slider->value();

    // Given "up" and "forward" directions, we can easily define a "left"
    // direction using their cross product.
    Vec3f left = cross(up, forward);

	bool update_dir = false;

    if( Fl::event_state(FL_ALT) )
        // With ALT held down, move camera parallel to image plane
        switch( key )
        {
        case FL_Left:   eye += move_step * left; break;
        case FL_Right:  eye -= move_step * left; break;

        case FL_Up:     eye += move_step * up; break;
        case FL_Down:   eye -= move_step * up; break;

        default: return false;
        }
    else
        switch( key )
        {
	 // Rotate viewing direction left & right
		case FL_Left: 
			theta += 5.0f;
			update_dir = true;
			break;
        case FL_Right:
			theta -= 5.0f;
			update_dir = true;
			break;

	// Move forward & backward
        case FL_Up:     eye += move_step * forward; break;
        case FL_Down:   eye -= move_step * forward; break;

	// Spin camera about the viewing axis
        //case FL_Page_Up:    up = rotate_direction(up, 10, forward); break;
        //case FL_Page_Down:  up = rotate_direction(up, -10, forward); break;

	// Rotate
		case '8': 
			phi -= 10.0f;
			update_dir = true;
			break;
		case '2': 
			phi += 10.0f;
			update_dir = true;
			break;

	// Debug subdivision
		case '.': 		
			terrain.subdiv_level = min(terrain.subdiv_level+1, terrain.max_subdiv_level);
			break;		
		case ',':
			terrain.subdiv_level = max(terrain.subdiv_level-1, 0); 			
			break;

	// Camera
		case ' ':
			update_terrain_viewpoint = !update_terrain_viewpoint;
			if (!update_terrain_viewpoint) {
				status("Not following camera viewpoint");
			}

        default: return false;
        }

	if (update_dir) {
		update_rotation();
		update_dir = false;
	}

	update_pos();	

    canvas->redraw();
    return true;
}

void GUI::cmdline_file(const char *filename)
{
    if( filename )
    {
	if( !terrain.z )
	{
	    cerr << "Reading terrain data from file: " << filename << endl;
	    terrain.z = read_image(filename);
	}
	else if( !terrain.texture )
	{
	    cerr << "Reading texture data from file: " << filename << endl;
	    terrain.texture = read_image(filename);
	}
    }
}

void GUI::add_upper_controls(int& yfill, const int pad) 
{
	int size_x = 280;
	int size_y = 16;

	// Add various sliders
	yfill += pad;
	tau_slider = new Fl_Value_Slider(20,yfill, size_x, size_y, "Error (pixels)");
	tau_slider->type(FL_HOR_NICE_SLIDER);
	tau_slider->minimum(1.0f);
	tau_slider->maximum(16.0f);	
	tau_slider->value(4.0f);

	tri_slider = new Fl_Value_Slider(20+size_x+20, yfill, size_x, size_y, "Triangle limit (thousands)");
	tri_slider->type(FL_HOR_NICE_SLIDER);
	tri_slider->minimum(1);
	tri_slider->maximum(256);
	tri_slider->value(32);

	int tau_yfill = tau_slider->h() + tau_slider->labelsize();
	int tri_yfill = tri_slider->h() + tri_slider->labelsize();
	yfill += max(tau_yfill, tri_yfill);

	height_slider = new Fl_Value_Slider(20, yfill, size_x, size_y, "Fly height");
	height_slider->type(FL_HOR_NICE_SLIDER);
	height_slider->minimum(1);
	height_slider->maximum(256);
	height_slider->value(16);

	speed_slider = new Fl_Value_Slider(20+size_x+20, yfill, size_x, size_y, "Speed (meters/frame)");
	speed_slider->type(FL_HOR_NICE_SLIDER);
	speed_slider->minimum(0.0f);
	speed_slider->maximum(64.0f);
	speed_slider->value( 1.0f );

	int height_yfill = height_slider->h() + height_slider->labelsize();
	int speed_yfill = speed_slider->h() + speed_slider->labelsize();
	yfill += max(height_yfill, speed_yfill);
}

const char* choose_file(const char* msg, const char* pattern, const char* default_fname)
{
	//return default_fname;
	return fl_file_chooser(msg, pattern, default_fname);
}

int approx_lg(int n)
{
	int lgn = 0;
	while (n > 0) {
		n >>= 1;
		lgn++;
	}

	return lgn;
}

void init_terrain()
{	
    // Dimensions of Grand Canyon data set
    gui.terrain.pixel_spacing = 1.0f;  // pixels are 60m apart
    gui.terrain.height_unit = 10.004f / 60.0f;  // each gray level is 10.004 meters
    gui.terrain.height_stretch = 3.6f;  // artificially stretch z (for looks)
	gui.terrain.height_factor = gui.terrain.height_unit * gui.terrain.height_stretch;

    int width = gui.terrain.z->width();
    int height = gui.terrain.z->height();	
	
	gui.terrain.inv_h = 1.0f / (float)(height - 1);
	gui.terrain.inv_w = 1.0f / (float)(width - 1);

	float downsample_factor = 4096.0f / (float)(width-1);
    gui.terrain.pixel_spacing *= downsample_factor;	

	gui.terrain.init();
}

void init_gui()
{	
	gui.toplevel->label("Terrifly");

	int width = gui.terrain.z->width();
    int height = gui.terrain.z->height();

    Vec3f origin = gui.terrain.pixel_to_space(0,0);
    Vec3f opposite = gui.terrain.pixel_to_space(width-1, height-1);
    Vec3f diag = opposite - origin;
    diag[2] = 0;
    
    gui.zfar = 2*unitize(diag);

    gui.move_step = gui.terrain.pixel_spacing;
	gui.speed_slider->value(gui.move_step);
	
    gui.eye = origin + Vec3f(0, 0, 100);

    gui.up = Vec3f(0, 0, 1);

	gui.update_rotation();

    gui.add_toggle_menu("&Terrain/Draw Texture", FL_CTRL+'x',
			gui.will_draw_texture);

    gui.add_toggle_menu("&Terrain/Draw Wireframe", FL_CTRL+'w',
			gui.will_draw_wireframe);

	gui.add_toggle_menu("&Terrain/View Frustum Culling", FL_CTRL+'c',
			gui.terrain.do_view_frustum_culling);

	gui.add_toggle_menu("&Terrain/Adaptive Tau", FL_CTRL+'p',
			adaptive_tau);

	gui.add_toggle_menu("&Terrain/Guarantee triangle limit", FL_CTRL+'l',
			guarantee_tri_limit);

	gui.add_toggle_menu("&Drawing/GPS", FL_CTRL+'p',
		gps);

	gui.add_toggle_menu("&Drawing/Compensate for Y-axis rotation", FL_CTRL+'r',
			rot_compensate);

	gui.add_toggle_menu("&Drawing/Compensate for X-axis rotation", FL_CTRL+'r',
			phi_compensate);

	gui.add_toggle_menu("&Drawing/Draw error lines", FL_CTRL+'e',
			display_error_lines);

	gui.add_toggle_menu("&Misc/Run timedemo", FL_CTRL+'t',
			timedemo);

	gui.add_toggle_menu("&Misc/Recalc Active", FL_CTRL+'r',
			recalc_active);

	gui.add_toggle_menu("&Misc/Compiled Vertex Arrays", FL_CTRL+'v',
			use_vertex_array);

	gui.add_toggle_menu("&Misc/Ground Walk", FL_CTRL+'g',
			ground_walk);

	gui.add_toggle_menu("&Misc/Fly", FL_CTRL+'l',
			fly);

	gui.add_toggle_menu("&Misc/Ground Collision", FL_CTRL+'g',
			ground_collide);
}

int main(int argc, char *argv[])
{	
	cout << "\nTerrify, by Leslie Wu\n" << endl;

    gui.initialize(argc, argv);

	// read terrain z data set
    if( !gui.terrain.z )
    {		
		const char *zname = choose_file("Select height data:", "z*.tif",			
			"z4096_gcanyon.tif");

		if( !zname ) exit(0);
		
		cout << "Reading height map " << zname << "... ";
		gui.terrain.z = read_image(zname);
		cout << "done" << endl;
		
		if (!gui.terrain.z) {
			error_popup("\tInvalid terrain image");
		}
    }

	int w = gui.terrain.z->width();
	int h = gui.terrain.z->height();
	int min_dim = min(w, h);	
	gui.terrain.max_subdiv_level = 2*(approx_lg(min_dim)-1);
	gui.terrain.subdiv_level = gui.terrain.max_subdiv_level;

	if (fl_ask("Precalculate terrain image data set again?")) {
		precalc_terrain(*gui.terrain.z);		
	}

	// read pre-calculated bounds
	gui.terrain.errorImg = read_image("errorImg.tif");
	if (!gui.terrain.errorImg) {
		error_popup("Couldn't read pre-calculated error bounds");
	}

	gui.terrain.radiusImg.init("radii.dat", w, h);

    if( !gui.terrain.texture )
    {
		const char *rgbname =
			choose_file("Select texture data:", "rgb*.tif",					
		            "rgb4096_gcanyon.tif");
		if (!rgbname) {
			error_popup("Texture file not chosen.");
		}
		cout << "Reading texture file " << rgbname << "... ";
		gui.terrain.texture = read_image(rgbname);
		cout << "done" << endl;
		if (!gui.terrain.texture) {
			error_popup("Invalid terrain texture");
		}
    }

	cout << "Initializing GUI... ";

	init_terrain();
	init_gui();

	cout << "done" << endl;

    return gui.run();
}
