#ifdef GL_ES
precision mediump float;
#endif

uniform float u_time;
uniform vec2 u_resolution;
uniform vec2 u_mouse;

#define num_objs 20
#define num_samples 10
#define num_bounces 5
#define check_num 0.01
float g_seed = 0.0;

//
// Base object types
//

struct ray{
	vec3 origin;
	vec3 dir;
};

struct camera{
	ray transform;
	float focal_length;
};

struct sphere{
	vec3 center;
	float radius;

	vec3 color;
};

struct disk{
	vec3 center;
	vec3 normal;
	float radius;

	vec3 color;
};

struct hittable{
	// enum of what type of object this is
	// 0: sphere
	// 1: disk
	// 2: point_light
	int type;

	vec3 center;
	vec3 normal;

	float radius;

	vec3 color;
};

struct hit{
	float dist;
	int type;

	vec3 pos;
	vec3 norm;
	vec3 ref;

	vec3 color;
};

//
// Utilities
//

// random for point
// https://stackoverflow.com/a/10625698
float random(vec2 pos){
    vec2 K1 = vec2(
        23.14069263277926, // e^pi (Gelfond's constant)
         2.665144142690225 // 2^sqrt(2) (Gelfondâ€“Schneider constant)
    );
    return fract( cos( dot(pos,K1) ) * 12345.6789 );
}

vec2 hash_vec2(vec2 pos){
    vec2 K1 = vec2(
        23.14069263277926, // e^pi (Gelfond's constant)
         2.665144142690225 // 2^sqrt(2) (Gelfondâ€“Schneider constant)
    );
	return vec2(
				fract( cos( dot(pos,K1) ) * 12345.6789 ),
				fract( cos( dot(pos,K1) ) * 34567.8912 )
			);
}

vec3 hash_vec3(vec2 pos){
    vec2 K1 = vec2(
        23.14069263277926, // e^pi (Gelfond's constant)
         2.665144142690225 // 2^sqrt(2) (Gelfondâ€“Schneider constant)
    );
	return vec3(
				fract( cos( dot(pos,K1) ) * 12345.6789 ),
				fract( cos( dot(pos,K1) ) * 34567.8912 ),
				fract( cos( dot(pos,K1) ) * 23456.7891 )
			);
}

vec3 random_in_unit_sphere(vec2 seed) {
	vec3 h = hash_vec3(seed) * vec3(2.,6.28318530718,1.)-vec3(1,0,0);
    float phi = h.y;
    float r = pow(h.z, 1./3.);
	return r * vec3(sqrt(1.-h.x*h.x)*vec2(sin(phi),cos(phi)),h.x);
}

vec3 random_on_hemisphere(vec3 normal, vec2 seed){
	vec3 base_vector = random_in_unit_sphere(seed);
	if(dot(base_vector, normal) > 0.0){
		return base_vector;
	} else{
		return -base_vector;
	}
}

vec4 color_mix(vec4 base_col, vec4 top_col, float top_factor){
	return (base_col * (1.0 - top_factor)) + (top_col * top_factor);
}

vec3 pos_along_ray(ray r, float t){
	return r.origin + (r.dir * t);
}

//
// Intersection Functions
//

float intersect_sphere(hittable obj, ray r){
	vec3 oc = obj.center - r.origin;
	float a = dot(r.dir, r.dir);
    float b = -2.0 * dot(r.dir, oc);
    float c = dot(oc, oc) - obj.radius*obj.radius;
    float discriminant = b*b - 4.0*a*c;

	// if we miss the sphere, return -1.0
	if(discriminant < 0.0){
		return -1.0;

	// else return where along the ray we hit the sphere
	} else{
		return (-b - sqrt(discriminant)) / (2.0 * a);
	}
}

float intersect_disk(hittable obj, ray r){
	float denom = dot(obj.normal, r.dir);

	// if the ray isn't parallel the plane
	if(abs(denom) > 0.00001){
		float t = dot(obj.center - r.origin, obj.normal) / denom;
		if(t < 0.0){
			return -1.0;
		} else{
			vec3 hit_point = pos_along_ray(r, t);
			float dist_from_disk_center = distance(hit_point, obj.center);
			
			if(dist_from_disk_center < obj.radius){
				return t;
			} else{
				return -1.0;
			}
		}
	}
}

hit intersect_world(ray r, hittable objs[num_objs]){

	bool hit_something = false;
	hit min_hit = hit(
				1000.0,			// dist (t)
				-1,				// type
				vec3(0.0),		// location
				vec3(0.0),		// normal
				vec3(0.0),		// reflect
				vec3(0.0)		// color
			);

	for(int obj = 0; obj < num_objs; obj++){
		
		float t = 0.0;
		if(objs[obj].type == 0 || objs[obj].type == 2){
			t = intersect_sphere(objs[obj], r);
		} else if(objs[obj].type == 1){
			t = intersect_disk(objs[obj], r);
		}
		if(t > 0.0001){
			hit_something = true;
			
			if(t < min_hit.dist){
				min_hit.dist = t;
				min_hit.pos = pos_along_ray(r, t);
				min_hit.color = objs[obj].color;
				min_hit.type = objs[obj].type;
				
				if(objs[obj].type == 0 || objs[obj].type == 2){
					min_hit.norm = min_hit.pos - objs[obj].center;
				} else if(objs[obj].type == 1){
					min_hit.norm = objs[obj].normal;
				}
				
				min_hit.ref = r.dir - 2.0*dot(r.dir, min_hit.norm)*min_hit.norm;
			}
		}
	}

	if(!hit_something){
		min_hit.dist = -1.0;;
	}
	
	return min_hit;
}


void mainImage(out vec4 frag_color, in vec2 frag_coord)
{
	// random seed
	g_seed = float(random(frag_coord))/float(0xffffffff)+u_time;
    
	// Normalized pixel coordinates (from 0 to 1)
    vec2 uv = frag_coord / u_resolution.xy;
	vec2 mouse_uv = u_mouse / u_resolution.xy;
	
	// Map the coordinates to have square pixels
	float smaller_side = min(u_resolution.x, u_resolution.y);
	vec2 mapped_coords = frag_coord / smaller_side;
	vec2 mapped_res = u_resolution / smaller_side;
	vec2 mouse_mapped = u_mouse / smaller_side;

	if(smaller_side == u_resolution.x){
		float larger_excess = (u_resolution.y/smaller_side) - 1.0;
		mapped_coords = vec2(mapped_coords.x, mapped_coords.y - (0.5 * larger_excess));
		mapped_res = vec2(mapped_res.x, mapped_res.y - (0.5 * larger_excess));
		mouse_mapped = vec2(mouse_mapped.x, mouse_mapped.y - (0.5 * larger_excess));
	} else if(smaller_side == u_resolution.y){
		float larger_excess = (u_resolution.x/smaller_side) - 1.0;
		mapped_coords = vec2(mapped_coords.x - (0.5 * larger_excess), mapped_coords.y);
		mapped_res = vec2(mapped_res.x - (0.5 * larger_excess), mapped_res.y);
		mouse_mapped = vec2(mouse_mapped.x - (0.5 * larger_excess), mouse_mapped.y);
	}

	// The base color of the output image
	vec4 base_color = vec4(0.0);

	// scene camera definition
	//vec3 cam_origin = vec3(0.0, 10.0, 40.0);
	vec3 cam_origin = vec3(0.0, 3.0, 40.0);
	camera cam = camera(
					ray(
						cam_origin,	// origin
						//vec3(0.0, 0.0, -1.0)	// direction
						-normalize(vec3(cam_origin.x, 0.0, cam_origin.z))
					),
					3.0							// focal length
				);
	//
	// direction to each pixel from the camera
	//

	// The direction to the current pixel
    vec3 pixel_dir = (
				(cam.transform.dir * cam.focal_length) + 
				vec3(2.0 * mapped_coords - 1.0, 0.0)
			);
	ray pixel_ray = ray(
				cam.transform.origin,
				normalize(pixel_dir)
			);

	// Background gradient
	// vec4 top_col = vec4(0.3, 0.3, 0.7, 1.0);
	// vec4 bot_col = vec4(0.8, 0.8, 1.0, 1.0);
	vec4 bot_col = vec4(vec3(0.05), 1.0);
	vec4 top_col = vec4(vec3(0.3), 1.0);
	vec4 gradient = color_mix(bot_col, top_col, 0.5 + (pixel_ray.dir.y * 0.5));

	// Define the scene
	
	hittable objs[num_objs];

	#define wall_size 50.0

	objs[0] = hittable(
				0,
				vec3(0.0, 1.0, 0.0),
				vec3(0.0, 0.0, 0.0),
				3.0,
				vec3(0.9, 0.9, 0.9)
			);
	objs[1] = hittable(							// mouse-tracker sphere
				2,
				vec3((mouse_mapped.x*2.0-1.0)*10.0,
					(mouse_mapped.y*2.0-1.0)*12.0,
					5.0),
				vec3(0.0, 0.0, 0.0),
				2.0,
				vec3(1.0, 1.0, 1.0)
			);
	objs[2] = hittable(					// floor
				1,						// type: disk
				vec3(0.0, -2.0, 0.0),	// center
				vec3(0.0, 1.0, 0.0),	// normal
				wall_size,					// radius
				vec3(0.9, 0.9, 0.9)		// color
			);
	objs[3] = hittable(					// left wall
				1,						// type: disk
				vec3(-20.0, 0.0, 0.0),	// center
				vec3(1.0, 0.0, 0.0),	// normal
				wall_size,					// radius
				vec3(1.0, 0.0, 0.0)		// color
			);
	// objs[4] = hittable(					// back wall
	// 			1,						// type: disk
	// 			vec3(0.0, 0.0, -20.0),	// center
	// 			vec3(0.0, 0.0, 1.0),	// normal
	// 			wall_size,					// radius
	// 			vec3(0.9, 0.9, 0.9)		// color
	// 		);
	objs[5] = hittable(					// right wall
				1,						// type: disk
				vec3(20.0, 0.0, 0.0),	// center
				vec3(-1.0, 0.0, 0.0),	// normal
				wall_size,					// radius
				vec3(0.0, 1.0, 0.0)		// color
			);
	objs[6] = hittable(					// ceiling
				1,						// type: disk
				vec3(0.0, 10.0, 0.0),	// center
				vec3(0.0, -1.0, 0.0),	// normal
				wall_size,					// radius
				vec3(0.9, 0.9, 0.9)		// color
			);

	for(int i = 7; i < num_objs; i++){
		vec3 rand_pos = hash_vec3(vec2(i+1, i+2)) * 30.0 - 15.0;
		vec3 rand_norm = hash_vec3(vec2(i*10, i*5)) * 45.0 - 22.0;

		int rand_type = int(random(vec2(i, i+1))*2.2);
		if(rand_type > 2){ rand_type == 2; }
		//int rand_type = 0;

		objs[i] = hittable(
					rand_type,
					vec3(rand_pos.x, rand_pos.y/10.0 + 2.0, rand_pos.z),
					vec3(rand_norm.x, rand_norm.y, rand_norm.z),
					random(vec2(i+1, i+2))*2.0+1.0,
					hash_vec3(vec2(i+3, i+4))
				);
	}

	//
    // Handle intersections
	//

	// handle bounces in a for loop because the current toolchain
	// doesn't support recursion in functions	
	vec3 total_color = vec3(0.0);

	for(int rand = 0; rand < num_samples; rand++){
		vec3 sample_color = vec3(0.0);
		ray cur_ray = pixel_ray;
		for(int bounce = 0; bounce < num_bounces; bounce++){
			hit ray_hit = intersect_world(cur_ray, objs);

			//
			// handle lighting/color
			//

			// if we hit an object, add the color
			if(ray_hit.dist > 0.0){
				if(ray_hit.type == 2){
					sample_color += ray_hit.color;
					break;
				} else{
					sample_color *= (0.001);
					sample_color += ((1.0/float(bounce+3)) * ray_hit.color);
					//sample_color += (0.1 * ray_hit.color);
				}

			// else add background color and break
			} else{
				break;
			}
			
			//
			// change cur_ray for the next bounce
			//
			
			// check if light could be found from a bounce
			vec3 light_dir = vec3(0.0);
			for(int obj = 0; obj < num_objs; obj++){
				if(objs[obj].type == 2){
					ray potential_ray = ray(ray_hit.pos, objs[obj].center - ray_hit.pos);
					hit light_hit = intersect_world(potential_ray, objs);

					if(light_hit.type == 2){
						light_dir = potential_ray.dir;
					}
				}
			}

			cur_ray = ray(
				ray_hit.pos,
				0.005*ray_hit.norm +
					random_on_hemisphere(
							ray_hit.norm,
							vec2(uv.x + float(rand)*g_seed, uv.y+float(rand+1)*g_seed)
						)
			);
		}
		
		if(sample_color != vec3(0.0)){
			total_color += sample_color;
		}
	}


	if(total_color == vec3(0.0)){
		total_color = gradient.xyz;
	} else{
		total_color /= float(num_samples);
	}
	
	frag_color = vec4(total_color, 1.0);
}

void main() {
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
