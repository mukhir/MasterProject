uniform vec3 cameravector;
uniform vec4 KFOV;
//uniform float minDis;
//uniform float disConstant;

//varying float depth;


void main()
{	
	// Shading
	
	vec3 N = normalize( gl_NormalMatrix * gl_Normal );
	
	vec3 P = ( gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0) ).xyz;
	
	vec3 V = normalize(-gl_Vertex.xyz);
	vec3 L = normalize(gl_LightSource[0].position.xyz);

	if( dot(V, N) < 0.0 )
		N = -N;
	//if( dot(N, L) < 0.0 )
	//	N = -N;

	vec3 H = normalize(L+V);
	
	float d = max( abs(dot(N, L)), 0.0 );
	float s = 0.0;
	//if( d > 0.0 )
	s = pow( max( abs(dot(N,H)), 0.0 ), gl_FrontMaterial.shininess );
	
	H = normalize(-L+V);
	s += pow( max( abs(dot(N,H)), 0.0 ), gl_FrontMaterial.shininess );
	
	vec4 color = gl_FrontMaterial.ambient * gl_LightSource[0].ambient + 
			gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse * d + 
			gl_FrontMaterial.specular * gl_LightSource[0].specular * s;

	gl_FrontColor = color;

	
	vec3 diff = gl_Vertex.xyz - cameravector;
	float dis = length( diff );
	float out_size;
	
	
	float att = 1.0 + 0.000001 * dis + 0.000000015 * dis * dis;	
	

	
	out_size =  max(1.0, (KFOV.y/abs(P.z)) * gl_Vertex.w );
	gl_PointSize = 1.1 * out_size/2.0;
		
	//Should be overwritten by point sprite
	//gl_TexCoord[0] = gl_MultiTexCoord0;
	
	//depth = abs((dis - minDis) * disConstant);
	//depth = abs((P.z - 0.1)/(1000.0 - 0.1));
	
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
}
