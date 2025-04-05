#ifdef GL_ES
precision mediump float;
#endif

uniform float u_time;
uniform vec2 u_resolution;
uniform vec2 u_mouse;

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

//
// Utilities
//

// random for point
// https://stackoverflow.com/a/10625698
float random(vec2 p){
    vec2 K1 = vec2(
        23.14069263277926, // e^pi (Gelfond's constant)
         2.665144142690225 // 2^sqrt(2) (Gelfondâ€“Schneider constant)
    );
    return fract( cos( dot(p,K1) ) * 12345.6789 );
}

float linear_congruential_generator(float prev){
	float modulus = 100.0;
	float multiplier = 13.0;
	float increment = 2.0;
	return mod(multiplier * prev + increment, modulus);
}

float rand_sequence(int n){
	float seed = 123.0;
	float cur = seed;
	int i = 0;
	for(int i = 0; i < 100; i++){
		cur = linear_congruential_generator(cur);
		if(i >= n){
			break;
		}
	}
	return cur;
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

float intersect_sphere(sphere obj, ray r){
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

float intersect_disk(disk obj, ray r){
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

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord / u_resolution.xy;
	vec2 mouse_uv = u_mouse / u_resolution.xy;
	
	// Map the coordinates to have square pixels
	float smaller_side = min(u_resolution.x, u_resolution.y);
	vec2 mapped_coords = fragCoord / smaller_side;
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
	
	// direction to each pixel from the camera
	
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
	vec4 top_col = vec4(0.3, 0.3, 0.7, 1.0);
	vec4 bot_col = vec4(0.8, 0.8, 1.0, 1.0);
	vec4 gradient = color_mix(bot_col, top_col, 0.5 + (pixel_ray.dir.y * 0.5));

	// Define the scene
	
	#define num_spheres 30
	sphere spheres[num_spheres];

	spheres[0] = sphere(
				vec3(0.0, -100.0, 0.0),
				100.0,
				vec3(0.9, 0.9, 0.9)
			);
	spheres[1] = sphere(
				vec3(0.0, 1.0, 0.0),
				3.0,
				vec3(0.9, 0.9, 0.9)
			);
	spheres[2] = sphere(							// mouse-tracker sphere
				vec3((mouse_mapped.x*2.0-1.0)*20.0,
					(mouse_mapped.y*2.0-1.0)*20.0,
					5.0),
				2.0,
				vec3(1.0, 0.0, 0.0)
			);
	// spheres[3] = sphere(							// light sphere
	// 			vec3(10.0, 10.0, 0.0),
	// 			3.0,
	// 			vec3(1.0, 1.0, 1.0)
	// 		);

	for(int i = 4; i < num_spheres; i++){
		vec3 rand_pos = hash_vec3(vec2(i, i))*45.0-22.0;
		spheres[i] = sphere(
					vec3(rand_pos.x, rand_pos.y/10.0 + 2.0, rand_pos.z),
					random(vec2(i+1, i+1))*2.0+1.0,
					hash_vec3(vec2(i+2, i+2))
				);
	}

    // Handle intersections
	//
	
	fragColor = gradient;

	// handle bounces in a for loop	
	ray cur_ray = pixel_ray;
	vec3 total_color = vec3(0.0, 0.0, 0.0);
	
	for(int bounce = 0; bounce < 2; bounce++){
		bool hit = false;
		for(int rand = 0; rand < 100; rand++){
			float min_t = 1000.0;
			vec3 min_hit_location = vec3(0.0);
			vec3 min_hit_normal = vec3(0.0);
			vec3 min_hit_color = gradient.xyz;
			vec3 min_hit_reflect = vec3(0.0);
			for(int i = 0; i < num_spheres; i++){
				// t: the position along a ray where it intersected
				// 		with an object
				float t = intersect_sphere(spheres[i], cur_ray);
				if(t > 0.0001){
					hit = true;
					if(t < min_t){
						min_t = t;
						min_hit_location = pos_along_ray(cur_ray, t);
						min_hit_normal = min_hit_location - spheres[i].center;
						min_hit_color = spheres[i].color;
						min_hit_reflect = cur_ray.dir - 2.0*dot(cur_ray.dir, min_hit_normal)*min_hit_normal;
					}
				}
			}
			cur_ray = ray(
						min_hit_location,
						0.005*min_hit_normal +
							random_on_hemisphere(
									min_hit_normal,
									vec2(uv.x + float(rand), uv.y+float(rand))
									//(fragCoord.x*fragCoord.y)+float(rand) // random seed
									// random(fragCoord)+float(rand)
								)
					);
			
			if(bounce == 0){
				if(min_hit_color.x+min_hit_color.y+min_hit_color.z >= 2.9){
					hit = false;
				}
				total_color = min_hit_color;
				break;
			} else{	
				total_color *= 0.7;
				total_color += (0.2 * min_hit_color);
				if(min_hit_color.x+min_hit_color.y+min_hit_color.z >= 2.9){
					total_color = min_hit_color*2.0;
				}
			}

			if(min_t == 1000.0){
				break;
			}
		}
		if(!hit){
			break;
		}

	}

	if(total_color != vec3(0.0)){
		fragColor = vec4(total_color, 1.0);
	}
	
	//float i = random(uv);
	//fragColor = vec4(vec3(i), 1.0);
}

void main() {
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
