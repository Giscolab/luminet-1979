
(async function(){
  const canvas = document.getElementById('c');
  const gl = canvas.getContext('webgl2');
  if(!gl){ alert('WebGL2 requis'); return; }

  // Vertex shader
  const vsSource = `#version 300 es
  in vec2 aPos;
  void main(){ gl_Position = vec4(aPos,0.0,1.0); }`;

  // Fragment shader
  const fsSource = `#version 300 es
  precision highp float;
  out vec4 fragColor;
  uniform vec2 iResolution;
  uniform float iTime;
  uniform float iIncl;
  uniform float iRout;
  uniform float iSpeed;
  uniform float iBloom;

  #define PI 3.141592653589793
  #define ISCO 6.0
  #define BC 5.196152

  float hash(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }
  float r_lens_primary(float b){ return b+2.0+4.0/b; }
  float r_lens_secondary(float b){ float d=abs(b-BC); return ISCO+60.0*exp(-d*1.8); }
  float redshift(float r,float phi,float incl){ float v=1.0/sqrt(max(r-3.0,0.001)); return max((1.0+v*cos(phi)*sin(incl))/sqrt(max(1.0-3.0/r,0.001)),0.001); }
  float emissivity(float r,float phi,float t){
    if(r<ISCO) return 0.0;
    float base=(1.0-sqrt(ISCO/r))*pow(r,-3.0);
    float spin=fract(phi/(2.0*PI)+t*0.1*iSpeed);
    float flare=pow(sin(spin*20.0+r*0.5)*0.5+0.5,4.0)*0.3;
    return base*(1.0+flare);
  }
  float dust(vec2 pos,float t){ float r=pos.x; float phi=pos.y; float rot=t*0.15*iSpeed*pow(r,-1.2); float angle=phi+rot; vec2 q=vec2(floor(r*0.8),floor(angle*12.0/PI)); float s=hash(q); float c=hash(vec2(floor(r*2.3),floor(angle*28.0))); return 0.55+0.28*s+0.17*c; }
  vec3 tempColor(float T){ T=clamp(T,0.0,1.0); vec3 c0=vec3(0.5,0.0,0.0), c1=vec3(1.0,0.2,0.0), c2=vec3(1.0,0.72,0.2), c3=vec3(1.0,0.98,0.9), c4=vec3(0.72,0.85,1.0); if(T<0.25) return mix(c0,c1,T/0.25); else if(T<0.55) return mix(c1,c2,(T-0.25)/0.30); else if(T<0.80) return mix(c2,c3,(T-0.55)/0.25); else return mix(c3,c4,(T-0.80)/0.20); }
  vec3 aces(vec3 x){ const float a=2.51,b=0.03,c=2.43,d=0.59,e=0.14; return clamp((x*(a*x+b))/(x*(c*x+d)+e),0.0,1.0); }
  vec3 bloom(vec3 col,float intensity){ float l=dot(col,vec3(0.2126,0.7152,0.0722)); if(l>0.8) col+=vec3(0.4,0.3,0.2)*(l-0.8)*2.0; return col; }

  void main(){
    vec2 uv=(gl_FragCoord.xy-0.5*iResolution.xy)/iResolution.y*30.0;
    float b=length(uv);
    float alpha=atan(uv.y,uv.x);

    float twink=sin(iTime*2.0+gl_FragCoord.x*0.01)*0.2+0.8;
    float stars=step(0.985,hash(floor(uv*18.0)))*(0.3+0.7*hash(uv*37.3))*twink;
    vec3 bg=vec3(stars*0.8,stars*0.85,stars);

    if(b<BC){
      float edge=smoothstep(BC,BC-0.3,b);
      bg=mix(bg*0.05,vec3(0.0),edge);
      float glow=smoothstep(BC,BC-0.5,b)*0.2;
      bg+=vec3(0.5,0.3,0.1)*glow*(0.8+0.4*sin(iTime*5.0));
      fragColor=vec4(bg,1.0); return;
    }

    vec3 col=bg;
    col+=vec3(0.9,0.6,0.3)*exp(-pow((b-BC)*1.5,2.0))*0.06;

    float r1=r_lens_primary(b);
    float z1=redshift(r1,alpha,iIncl);
    float flux1=emissivity(r1,alpha,iTime)*pow(1.0/z1,4.0);
    flux1*=dust(vec2(r1,alpha),iTime)*200.0;
    vec3 c1=tempColor(clamp(0.5+(1.0/z1-1.0)*1.2,0.0,1.0))*flux1;

    vec3 c2=vec3(0.0);
    if(b<25.0){
      float r2=r_lens_secondary(b);
      float z2=redshift(r2,alpha+PI,iIncl);
      float flux2=emissivity(r2,alpha+PI,iTime)*pow(1.0/z2,4.0);
      flux2*=dust(vec2(r2,alpha+PI),iTime)*60.0;
      c2=tempColor(clamp(0.5+(1.0/z2-1.0)*1.2,0.0,1.0))*flux2*0.3;
    }
    col+=c1+c2;
    col=aces(col*1.2);
    if(iBloom>0.5) col=bloom(col,1.0);
    float vignette=1.0-0.2*length(uv/30.0);
    col*=vignette;
    col=pow(col,vec3(1.0/2.2));
    fragColor=vec4(col,1.0);
  }`;

  function compileShader(type, source){
    const shader=gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS)){
      console.error(gl.getShaderInfoLog(shader));
    }
    return shader;
  }

  const vs=compileShader(gl.VERTEX_SHADER, vsSource);
  const fs=compileShader(gl.FRAGMENT_SHADER, fsSource);
  const prog=gl.createProgram();
  gl.attachShader(prog, vs);
  gl.attachShader(prog, fs);
  gl.linkProgram(prog);
  gl.useProgram(prog);

  const vao=gl.createVertexArray();
  gl.bindVertexArray(vao);
  const buf=gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,1,1]),gl.STATIC_DRAW);
  const posLoc=gl.getAttribLocation(prog,'aPos');
  gl.enableVertexAttribArray(posLoc);
  gl.vertexAttribPointer(posLoc,2,gl.FLOAT,false,0,0);

  const uRes=gl.getUniformLocation(prog,'iResolution');
  const uTime=gl.getUniformLocation(prog,'iTime');
  const uIncl=gl.getUniformLocation(prog,'iIncl');
  const uRout=gl.getUniformLocation(prog,'iRout');
  const uSpeed=gl.getUniformLocation(prog,'iSpeed');
  const uBloom=gl.getUniformLocation(prog,'iBloom');

  let incl=83*Math.PI/180, rout=36.0, speed=1.0, bloomMode=0;

  function resize(){
    const dpr=window.devicePixelRatio||1;
    canvas.width=Math.floor(canvas.clientWidth*dpr);
    canvas.height=Math.floor(canvas.clientHeight*dpr);
    gl.viewport(0,0,canvas.width,canvas.height);
  }
  window.addEventListener('resize',resize);
  resize();

  let lastTime=0, frames=0, fpsTimer=0;
  const fpsSpan=document.getElementById('fps');

  function render(now){
    now*=0.001;
    frames++;
    if(now-fpsTimer>=1.0){
      fpsSpan.textContent=Math.round(frames/(now-fpsTimer))+' fps';
      frames=0;
      fpsTimer=now;
    }
    gl.uniform2f(uRes,canvas.width,canvas.height);
    gl.uniform1f(uTime,now);
    gl.uniform1f(uIncl,incl);
    gl.uniform1f(uRout,rout);
    gl.uniform1f(uSpeed,speed);
    gl.uniform1f(uBloom,bloomMode);
    gl.drawArrays(gl.TRIANGLE_STRIP,0,4);
    requestAnimationFrame(render);
  }
  requestAnimationFrame(render);
})();

