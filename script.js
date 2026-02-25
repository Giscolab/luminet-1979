(async function () {
  const canvas = document.getElementById('c');
  const gl = canvas.getContext('webgl2');
  if (!gl) {
    alert('WebGL2 requis');
    return;
  }

  const vsSource = `#version 300 es
  in vec2 aPos;
  void main(){ gl_Position = vec4(aPos,0.0,1.0); }`;

  const fsSource = `#version 300 es
  precision highp float;
  out vec4 fragColor;
  uniform vec2 iResolution;
  uniform float iTime;
  uniform float iIncl;
  uniform float iRout;
  uniform float iSpeed;
  uniform float iBloom;
  uniform vec3 iCamPos;
  uniform vec3 iCamRight;
  uniform vec3 iCamUp;
  uniform vec3 iCamForward;
  uniform float iDiskHalfThickness;
  uniform float iBend;
  uniform float iSpin;
  uniform float iQuality;
  uniform float iFovX;
  uniform vec3 iCamVel;
  uniform float iFogStrength;
  uniform float iDiskSense;
  uniform float iEmissionMode;
  uniform float iMaxDiskCrossings;
  uniform float iLodAuto;
  uniform float iLodMin;
  uniform float iLodMax;
  uniform float iCamSpeedFiltered;

  #define PI 3.141592653589793
  #define BC 5.196152
  #define MAX_STEPS 300
  #define SKY_RADIUS 180.0

  float hash(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }
  float iscoRadius(float a, bool prograde){
    float aa = clamp(a, -0.999, 0.999);
    float z1 = 1.0 + pow(1.0 - aa * aa, 1.0 / 3.0) * (pow(1.0 + aa, 1.0 / 3.0) + pow(1.0 - aa, 1.0 / 3.0));
    float z2 = sqrt(max(3.0 * aa * aa + z1 * z1, 1e-5));
    float sgn = prograde ? 1.0 : -1.0;
    float inside = max((3.0 - z1) * (3.0 + z1 + 2.0 * z2), 0.0);
    return 3.0 + z2 - sgn * sqrt(inside);
  }
  float horizonRadius(float a){
    float aa = clamp(a, -0.999, 0.999);
    return 1.0 + sqrt(max(1.0 - aa * aa, 0.0));
  }
  float redshift(float r,float phi,float incl){ float v=1.0/sqrt(max(r-3.0,0.001)); return max((1.0+v*cos(phi)*sin(incl))/sqrt(max(1.0-3.0/r,0.001)),0.001); }
  float legacyEmissivity(float r,float phi,float t,float rISCO){
    if(r<rISCO || r>iRout) return 0.0;
    float base=(1.0-sqrt(rISCO/r))*pow(r,-3.0);
    float spin=fract(phi/(2.0*PI)+t*0.1*iSpeed);
    float flare=pow(sin(spin*20.0+r*0.5)*0.5+0.5,4.0)*0.3;
    return base*(1.0+flare);
  }

  float ntFlux(float r, float a, float rISCO){
    if(r<rISCO || r>iRout) return 0.0;
    float rr = max(r, rISCO + 1e-3);
    float rin = max(rISCO, 1.0);
    float edge = max(1.0 - sqrt(rin / rr), 0.0);
    float radial = pow(rr, -3.0);
    float spinBoost = 1.0 + 0.2 * a / max(pow(rr, 1.5), 1.0);
    float outerTaper = 1.0 / (1.0 + pow(rr / max(iRout * 1.12, rin + 0.5), 4.0));
    float flux = edge * radial * max(spinBoost, 0.15) * outerTaper;
    float floorFlux = 3e-4 * pow(rin, -3.0);
    return max(flux, floorFlux);
  }

  float emissivity(float r,float phi,float t,float rISCO){
    if(iEmissionMode < 0.5) return legacyEmissivity(r, phi, t, rISCO);
    float phiWarp = sin(phi * 2.0 + t * 0.25 * iSpeed) * 0.02;
    float warpedR = max(r + phiWarp * r, rISCO + 1e-3);
    return ntFlux(warpedR, iSpin, rISCO);
  }
  float dust(vec2 pos,float t){
    float r=pos.x; float phi=pos.y;
    float rot=t*0.15*iSpeed*pow(max(r,0.001),-1.2);
    float angle=phi+rot;
    vec2 q=vec2(floor(r*0.8),floor(angle*12.0/PI));
    float s=hash(q);
    float c=hash(vec2(floor(r*2.3),floor(angle*28.0)));
    return 0.55+0.28*s+0.17*c;
  }
  vec3 tempColor(float T){ T=clamp(T,0.0,1.0); vec3 c0=vec3(0.5,0.0,0.0), c1=vec3(1.0,0.2,0.0), c2=vec3(1.0,0.72,0.2), c3=vec3(1.0,0.98,0.9), c4=vec3(0.72,0.85,1.0); if(T<0.25) return mix(c0,c1,T/0.25); else if(T<0.55) return mix(c1,c2,(T-0.25)/0.30); else if(T<0.80) return mix(c2,c3,(T-0.55)/0.25); else return mix(c3,c4,(T-0.80)/0.20); }
  vec3 blackbodyApprox(float tNorm){
    float t = clamp(tNorm, 0.0, 1.35);
    vec3 red = vec3(1.0, 0.34, 0.08);
    vec3 amber = vec3(1.0, 0.69, 0.32);
    vec3 warmWhite = vec3(1.0, 0.94, 0.82);
    vec3 blueWhite = vec3(0.74, 0.85, 1.0);
    if(t < 0.45) return mix(red, amber, t / 0.45);
    if(t < 0.85) return mix(amber, warmWhite, (t - 0.45) / 0.40);
    return mix(warmWhite, blueWhite, (t - 0.85) / 0.50);
  }

  vec3 diskColor(float flux, float z, float dustTerm, float rISCO){
    if(iEmissionMode < 0.5){
      float legacyTemp = clamp(0.5 + (1.0 / max(z, 0.1) - 1.0) * 1.2, 0.0, 1.0);
      return tempColor(legacyTemp) * flux * dustTerm;
    }
    float refFlux = max(ntFlux(rISCO * 1.5, iSpin, rISCO), 1e-6);
    float Tcol = pow(max(flux, 1e-8) / refFlux, 0.25);
    float observedT = Tcol / max(z, 0.32);
    vec3 bb = blackbodyApprox(observedT);
    return bb * flux * dustTerm;
  }
  vec3 aces(vec3 x){ const float a=2.51,b=0.03,c=2.43,d=0.59,e=0.14; return clamp((x*(a*x+b))/(x*(c*x+d)+e),0.0,1.0); }
  vec3 bloom(vec3 col,float intensity){ float l=dot(col,vec3(0.2126,0.7152,0.0722)); if(l>0.8) col+=vec3(0.4,0.3,0.2)*(l-0.8)*2.0*intensity; return col; }

  vec3 aberrateDirection(vec3 dir, vec3 vel){
    float beta = clamp(length(vel), 0.0, 0.92);
    if (beta < 1e-4) return normalize(dir);
    vec3 n = vel / beta;
    float gamma = inversesqrt(max(1.0 - beta * beta, 1e-4));
    float dPar = dot(dir, n);
    vec3 par = n * dPar;
    vec3 perp = dir - par;
    float denom = 1.0 + beta * dPar;
    return normalize((perp / gamma + par + beta * n) / max(denom, 1e-4));
  }

  float dopplerFactor(vec3 lightDir, vec3 vel){
    float beta = clamp(length(vel), 0.0, 0.92);
    if (beta < 1e-4) return 1.0;
    vec3 n = vel / beta;
    float gamma = inversesqrt(max(1.0 - beta * beta, 1e-4));
    float mu = dot(-normalize(lightDir), n);
    return 1.0 / max(gamma * (1.0 - beta * mu), 1e-3);
  }

  vec3 applyDopplerTint(vec3 col, float doppler){
    float blueShift = clamp(doppler - 1.0, 0.0, 1.2);
    float redShift = clamp(1.0 - doppler, 0.0, 1.2);
    vec3 blue = vec3(0.78, 0.9, 1.15);
    vec3 red = vec3(1.2, 0.74, 0.55);
    col *= mix(vec3(1.0), blue, blueShift * 0.6);
    col *= mix(vec3(1.0), red, redShift * 0.55);
    col *= mix(0.65, 1.55, clamp(doppler, 0.45, 1.55) - 0.45);
    return col;
  }

  vec3 geoAccel(vec3 pos, vec3 dir, float L2){
    float r = max(length(pos), 0.001);
    float grav = -1.5 * iBend * L2 / pow(r, 5.0);
    vec3 acc = grav * pos;

    float frameDrag = 0.75 * iSpin * L2 / pow(r, 4.0);
    vec3 omega = vec3(0.0, frameDrag, 0.0);
    acc += cross(omega, dir);
    return acc;
  }

  void rk4Step(inout vec3 pos, inout vec3 dir, float h, float L2){
    vec3 k1p = dir;
    vec3 k1v = geoAccel(pos, dir, L2);

    vec3 p2 = pos + 0.5 * h * k1p;
    vec3 v2 = normalize(dir + 0.5 * h * k1v);
    vec3 k2p = v2;
    vec3 k2v = geoAccel(p2, v2, L2);

    vec3 p3 = pos + 0.5 * h * k2p;
    vec3 v3 = normalize(dir + 0.5 * h * k2v);
    vec3 k3p = v3;
    vec3 k3v = geoAccel(p3, v3, L2);

    vec3 p4 = pos + h * k3p;
    vec3 v4 = normalize(dir + h * k3v);
    vec3 k4p = v4;
    vec3 k4v = geoAccel(p4, v4, L2);

    pos += h * (k1p + 2.0*k2p + 2.0*k3p + k4p) / 6.0;
    dir = normalize(dir + h * (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0);
  }

  vec3 skyColor(vec3 dir, float starLod, float jitterAmp){
    vec3 d = normalize(dir);
    float phi = atan(d.z, d.x);
    float theta = acos(clamp(d.y, -1.0, 1.0));

    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);
    vec2 jitter = vec2(
      hash(gl_FragCoord.xy + vec2(iTime * 61.7, iTime * 23.3)) - 0.5,
      hash(gl_FragCoord.yx + vec2(iTime * 17.1, iTime * 43.9)) - 0.5
    ) * jitterAmp;
    vec2 lodGrid = mix(vec2(120.0, 60.0), vec2(220.0, 110.0), starLod);
    vec2 g = floor((uv + jitter / lodGrid) * lodGrid);

    float cluster = hash(g);
    float stars = step(0.994, cluster) * (0.25 + 0.75 * hash(g * 0.73 + 9.1));

    float milkyBand = exp(-pow((uv.y - 0.5) * 7.0, 2.0));
    float dustBand = hash(g * vec2(1.3, 2.1) + 17.0) * milkyBand;

    vec3 nebula = mix(vec3(0.03, 0.035, 0.055), vec3(0.08, 0.06, 0.04), dustBand * 0.8);
    vec3 starTint = mix(vec3(0.85, 0.9, 1.0), vec3(1.0, 0.88, 0.72), hash(g + 3.7));

    return nebula + starTint * stars * (0.7 + 0.3 * milkyBand);
  }

  bool traceScene(vec3 ro, vec3 rd, float rISCO, float traceQuality, out vec3 skyDir, out float bMin, out bool swallowed, out float fogOptical, out float fogGlow, out vec3 diskAccum){
    vec3 pos = ro;
    vec3 dir = rd;
    float L2 = dot(cross(ro, rd), cross(ro, rd));
    float prevY = pos.y;
    float diskTransmittance = 1.0;
    int crossings = 0;
    int crossingBudget = int(clamp(floor(mix(1.0, iMaxDiskCrossings, traceQuality) + 0.5), 1.0, 8.0));
    bMin = 1e9;
    swallowed = false;
    skyDir = rd;
    fogOptical = 0.0;
    fogGlow = 0.0;
    diskAccum = vec3(0.0);

    for (int i = 0; i < MAX_STEPS; i++) {
      if (float(i) > mix(80.0, float(MAX_STEPS), traceQuality)) break;

      float r = length(pos);
      bMin = min(bMin, length(cross(pos, dir)));
      float rH = horizonRadius(iSpin);
      float asymH = rH * (1.0 - 0.12 * iSpin * clamp(pos.x / max(r, 1e-3), -1.0, 1.0));
      if (r < clamp(asymH, max(0.8 * rH, 1.05), 2.45)) {
        swallowed = true;
        return false;
      }
      if (r > SKY_RADIUS) {
        skyDir = normalize(pos);
        return false;
      }

      float stepFar = mix(0.45, 1.4, clamp((r - 8.0) / 55.0, 0.0, 1.0));
      float h = stepFar / mix(1.65, 0.68, traceQuality);
      h *= mix(1.0, 0.72, clamp(10.0 / max(r, 0.1), 0.0, 1.0));

      float diskBand = exp(-abs(pos.y) / max(iDiskHalfThickness * 1.8, 0.05));
      float rXZ = length(pos.xz);
      float radial = smoothstep(rISCO * 0.7, rISCO * 1.2, rXZ)
                   * (1.0 - smoothstep(iRout * 0.9, iRout * 1.4, rXZ));
      float hotInner = exp(-pow((length(pos.xz) - (rISCO + 1.0)) * 0.45, 2.0));
      float fogDensity = diskBand * radial;
      float fogMask = step(0.01, fogDensity);
      fogOptical += fogDensity * h * 0.13 * iFogStrength * fogMask;
      fogGlow += fogDensity * (0.16 + 0.95 * hotInner) * h * iFogStrength * fogMask;

      vec3 prevPos = pos;
      vec3 prevDir = dir;
      rk4Step(pos, dir, h, L2);

      if (sign(prevY) != sign(pos.y)) {
        float t = prevY / (prevY - pos.y);
        vec3 crossPos = mix(prevPos, pos, t);
        float rCross = length(crossPos.xz);
        if (abs(crossPos.y) <= iDiskHalfThickness && rCross > rISCO && rCross < iRout) {
          vec3 crossDir = normalize(mix(prevDir, dir, t));
          float alpha = atan(crossPos.z, crossPos.x);

          float z1 = redshift(rCross, alpha, iIncl);
          float F1 = emissivity(rCross, alpha, iTime, rISCO);
          float flux1 = F1 * pow(1.0 / z1, 4.0);
          float dust1 = dust(vec2(rCross, alpha), iTime) * mix(170.0, 220.0, iEmissionMode);
          vec3 c1 = diskColor(flux1, z1, dust1, rISCO);

          float ghostOffset = 10.0 * exp(-abs(bMin - BC) * 1.7);
          float rGhost = clamp(rCross + ghostOffset, rISCO, iRout * 0.95);
          float z2 = redshift(rGhost, alpha + PI, iIncl);
          float F2 = emissivity(rGhost, alpha + PI, iTime, rISCO);
          float flux2 = F2 * pow(1.0 / z2, 4.0);
          float dust2 = dust(vec2(rGhost, alpha + PI), iTime) * mix(40.0, 55.0, iEmissionMode);
          float ghostWeight = 0.28 * smoothstep(0.0, 3.0, ghostOffset);
          vec3 c2 = diskColor(flux2, z2, dust2, rISCO) * ghostWeight;

          float thicknessFade = smoothstep(iDiskHalfThickness, 0.0, abs(crossPos.y));
          float dopplerAniso = 0.75 + 0.25 * clamp(dot(normalize(vec3(-crossPos.z,0.0,crossPos.x)), -crossDir), -1.0, 1.0);
          vec3 crossingCol = (c1 + c2) * thicknessFade * dopplerAniso;
          crossingCol = applyDopplerTint(crossingCol, dopplerFactor(crossDir, iCamVel));

          float radialNorm = clamp((rCross - rISCO) / max(iRout - rISCO, 1e-3), 0.0, 1.0);
          float localTau = mix(1.2, 0.18, radialNorm) * clamp(0.35 / max(abs(crossDir.y), 0.12), 0.7, 2.9);
          localTau *= mix(0.85, 1.25, iEmissionMode);
          float absorbWeight = 1.0 - exp(-localTau);
          diskAccum += crossingCol * (diskTransmittance * absorbWeight);
          diskTransmittance *= exp(-localTau);

          crossings++;
          if (crossings >= crossingBudget || diskTransmittance < 0.03) {
            skyDir = normalize(pos);
            return crossings > 0;
          }
        }
      }
      prevY = pos.y;
    }
    skyDir = normalize(pos);
    return crossings > 0;
  }

  void main(){
    vec2 p=(gl_FragCoord.xy-0.5*iResolution.xy)/iResolution.y;
    float focal = 0.5 * (iResolution.x / iResolution.y) / tan(radians(iFovX) * 0.5);
    vec3 rdCam = normalize(vec3(p.x, p.y, focal));
    vec3 rayDir = normalize(iCamRight * rdCam.x + iCamUp * rdCam.y + iCamForward * rdCam.z);
    rayDir = aberrateDirection(rayDir, iCamVel);

    bool prograde = iDiskSense > 0.5;
    float rISCO = iscoRadius(iSpin, prograde);

    float speed = length(iCamVel);
    float lodSpeed = mix(speed, iCamSpeedFiltered, 0.82);
    float lodNorm = smoothstep(0.04, 1.10, lodSpeed);
    float autoQuality = mix(iLodMax, iLodMin, lodNorm);
    float traceQuality = mix(iQuality, autoQuality, iLodAuto);
    float starLod = clamp(traceQuality * 0.95 + 0.05, 0.0, 1.0);
    float jitterAmp = mix(0.95, 0.18, traceQuality);

    vec3 skyDir = vec3(0.0);
    vec3 diskAccum = vec3(0.0);
    float bMin;
    float fogOptical = 0.0;
    float fogGlow = 0.0;
    bool swallowed = false;
    bool seenDisk = traceScene(iCamPos, rayDir, rISCO, traceQuality, skyDir, bMin, swallowed, fogOptical, fogGlow, diskAccum);

    vec3 bg = skyColor(skyDir, starLod, jitterAmp);
    float skyD = dopplerFactor(skyDir, iCamVel);
    bg = applyDopplerTint(bg, skyD);

    if (swallowed) {
      float photonRing = exp(-pow((bMin - BC) * 3.5, 2.0));
      bg = vec3(0.0);
      bg += vec3(0.95, 0.66, 0.35) * photonRing * 0.45;
      bg = applyDopplerTint(bg, dopplerFactor(rayDir, iCamVel));
      fragColor = vec4(pow(bg, vec3(1.0 / 2.2)), 1.0);
      return;
    }

    if (!seenDisk && bMin < BC) {
      float edge=smoothstep(BC,BC-0.5,bMin);
      bg=mix(bg*0.03,vec3(0.0),edge);
      float glow=smoothstep(BC,BC-0.8,bMin)*0.3;
      bg+=vec3(0.55,0.32,0.1)*glow*(0.8+0.4*sin(iTime*5.0));
      fragColor=vec4(bg,1.0); return;
    }

    vec3 col=bg;
    col+=vec3(0.9,0.6,0.3)*exp(-pow((bMin-BC)*1.3,2.0))*0.08;

    float fogTransmittance = exp(-fogOptical);
    col *= fogTransmittance;
    col += vec3(1.05, 0.72, 0.45) * fogGlow * 0.11;

    if(seenDisk) col += diskAccum;

    col=aces(col*1.2);
    if(iBloom>0.5) col=bloom(col,1.0);
    float vignette=1.0-0.2*length(p*1.1);
    col*=vignette;
    col=pow(col,vec3(1.0/2.2));
    fragColor=vec4(col,1.0);
  }`;

  function compileShader(type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      console.error(gl.getShaderInfoLog(shader));
    }
    return shader;
  }

  const vs = compileShader(gl.VERTEX_SHADER, vsSource);
  const fs = compileShader(gl.FRAGMENT_SHADER, fsSource);
  const prog = gl.createProgram();
  gl.attachShader(prog, vs);
  gl.attachShader(prog, fs);
  gl.linkProgram(prog);
  gl.useProgram(prog);

  const vao = gl.createVertexArray();
  gl.bindVertexArray(vao);
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
  const posLoc = gl.getAttribLocation(prog, 'aPos');
  gl.enableVertexAttribArray(posLoc);
  gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

  const uRes = gl.getUniformLocation(prog, 'iResolution');
  const uTime = gl.getUniformLocation(prog, 'iTime');
  const uIncl = gl.getUniformLocation(prog, 'iIncl');
  const uRout = gl.getUniformLocation(prog, 'iRout');
  const uSpeed = gl.getUniformLocation(prog, 'iSpeed');
  const uBloom = gl.getUniformLocation(prog, 'iBloom');
  const uCamPos = gl.getUniformLocation(prog, 'iCamPos');
  const uCamRight = gl.getUniformLocation(prog, 'iCamRight');
  const uCamUp = gl.getUniformLocation(prog, 'iCamUp');
  const uCamForward = gl.getUniformLocation(prog, 'iCamForward');
  const uDiskHalfThickness = gl.getUniformLocation(prog, 'iDiskHalfThickness');
  const uBend = gl.getUniformLocation(prog, 'iBend');
  const uSpin = gl.getUniformLocation(prog, 'iSpin');
  const uQuality = gl.getUniformLocation(prog, 'iQuality');
  const uFovX = gl.getUniformLocation(prog, 'iFovX');
  const uCamVel = gl.getUniformLocation(prog, 'iCamVel');
  const uFogStrength = gl.getUniformLocation(prog, 'iFogStrength');
  const uDiskSense = gl.getUniformLocation(prog, 'iDiskSense');
  const uEmissionMode = gl.getUniformLocation(prog, 'iEmissionMode');
  const uMaxDiskCrossings = gl.getUniformLocation(prog, 'iMaxDiskCrossings');
  const uLodAuto = gl.getUniformLocation(prog, 'iLodAuto');
  const uLodMin = gl.getUniformLocation(prog, 'iLodMin');
  const uLodMax = gl.getUniformLocation(prog, 'iLodMax');
  const uCamSpeedFiltered = gl.getUniformLocation(prog, 'iCamSpeedFiltered');

  const incVal = document.getElementById('incVal');
  const routVal = document.getElementById('routVal');
  const modeBadge = document.getElementById('mode-badge');
  const incSlider = document.getElementById('incSlider');
  const incSliderVal = document.getElementById('incSliderVal');
  const routSlider = document.getElementById('routSlider');
  const routSliderVal = document.getElementById('routSliderVal');
  const rotSlider = document.getElementById('rotSlider');
  const rotSliderVal = document.getElementById('rotSliderVal');
  const yawSlider = document.getElementById('yawSlider');
  const yawSliderVal = document.getElementById('yawSliderVal');
  const pitchSlider = document.getElementById('pitchSlider');
  const pitchSliderVal = document.getElementById('pitchSliderVal');
  const distSlider = document.getElementById('distSlider');
  const distSliderVal = document.getElementById('distSliderVal');
  const thickSlider = document.getElementById('thickSlider');
  const thickSliderVal = document.getElementById('thickSliderVal');
  const bendSlider = document.getElementById('bendSlider');
  const bendSliderVal = document.getElementById('bendSliderVal');
  const spinSlider = document.getElementById('spinSlider');
  const spinSliderVal = document.getElementById('spinSliderVal');
  const progradeToggle = document.getElementById('progradeToggle');
  const progradeVal = document.getElementById('progradeVal');
  const qualitySlider = document.getElementById('qualitySlider');
  const qualitySliderVal = document.getElementById('qualitySliderVal');
  const autoLodToggle = document.getElementById('autoLodToggle');
  const autoLodVal = document.getElementById('autoLodVal');
  const lodMinSlider = document.getElementById('lodMinSlider');
  const lodMinVal = document.getElementById('lodMinVal');
  const lodMaxSlider = document.getElementById('lodMaxSlider');
  const lodMaxVal = document.getElementById('lodMaxVal');
  const maxDiskCrossingsSlider = document.getElementById('maxDiskCrossingsSlider');
  const maxDiskCrossingsVal = document.getElementById('maxDiskCrossingsVal');
  const followOrbitToggle = document.getElementById('followOrbitToggle');
  const emissionModeToggle = document.getElementById('emissionModeToggle');
  const emissionModeVal = document.getElementById('emissionModeVal');
  const modeNormal = document.getElementById('modeNormal');
  const modeBloom = document.getElementById('modeBloom');
  const crosshair = document.getElementById('crosshair');

  let incl = 83 * Math.PI / 180;
  let rout = 36.0;
  let speed = 1.0;
  let bloomMode = 0;
  let camYaw = 0.0;
  let camPitch = 22.0;
  let camDist = 38.0;
  let diskHalfThickness = 0.9;
  let bend = 1.0;
  let spin = 0.2;
  let diskPrograde = true;
  let quality = 0.72;
  let autoLodEnabled = true;
  let lodMin = 0.40;
  let lodMax = 0.90;
  let filteredCamSpeed = 0.0;
  let maxDiskCrossings = 4;
  let ntEmissionMode = true;
  const fovX = 90.0;
  const fogStrength = 1.0;

  let followOrbit = true;

  function updateReadouts() {
    const inclDeg = (incl * 180 / Math.PI).toFixed(1);
    incVal.textContent = `${inclDeg}°`;
    routVal.textContent = rout.toFixed(1);
    incSlider.value = inclDeg;
    incSliderVal.textContent = `${inclDeg}°`;
    routSlider.value = rout.toFixed(1);
    routSliderVal.textContent = rout.toFixed(1);
    rotSlider.value = speed.toFixed(2);
    rotSliderVal.textContent = `${speed.toFixed(2)}x`;
    yawSlider.value = camYaw.toFixed(1);
    yawSliderVal.textContent = `${camYaw.toFixed(1)}°`;
    pitchSlider.value = camPitch.toFixed(1);
    pitchSliderVal.textContent = `${camPitch.toFixed(1)}°`;
    distSlider.value = camDist.toFixed(1);
    distSliderVal.textContent = camDist.toFixed(1);
    thickSlider.value = diskHalfThickness.toFixed(2);
    thickSliderVal.textContent = diskHalfThickness.toFixed(2);
    bendSlider.value = bend.toFixed(2);
    bendSliderVal.textContent = `${bend.toFixed(2)}x`;
    spinSlider.value = spin.toFixed(2);
    spinSliderVal.textContent = spin.toFixed(2);
    progradeToggle.checked = diskPrograde;
    progradeVal.textContent = diskPrograde ? "prograde" : "retrograde";
    qualitySlider.value = quality.toFixed(2);
    qualitySliderVal.textContent = quality.toFixed(2);
    autoLodToggle.checked = autoLodEnabled;
    autoLodVal.textContent = autoLodEnabled ? 'activé' : 'désactivé';
    lodMinSlider.value = lodMin.toFixed(2);
    lodMinVal.textContent = lodMin.toFixed(2);
    lodMaxSlider.value = lodMax.toFixed(2);
    lodMaxVal.textContent = lodMax.toFixed(2);
    maxDiskCrossingsSlider.value = String(maxDiskCrossings);
    maxDiskCrossingsVal.textContent = String(maxDiskCrossings);
    emissionModeToggle.checked = ntEmissionMode;
    emissionModeVal.textContent = ntEmissionMode ? "NT" : "stylized/legacy";
  }

  function setBloomMode(enabled) {
    bloomMode = enabled ? 1 : 0;
    modeNormal.classList.toggle('active', !enabled);
    modeBloom.classList.toggle('active', enabled);
  }

  function setModeBadge(text, cssClass = '') {
    modeBadge.textContent = text;
    modeBadge.className = cssClass;
  }

  function clamp(v, min, max) {
    return Math.max(min, Math.min(max, v));
  }

  function normalize(v) {
    const len = Math.hypot(v[0], v[1], v[2]) || 1;
    return [v[0] / len, v[1] / len, v[2] / len];
  }

  function cross(a, b) {
    return [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
  }

  function getSkimmingOrbitBasis(t) {
    const semiMajor = 12.0;
    const eccentricity = 0.68;
    const p = semiMajor * (1.0 - eccentricity * eccentricity);
    const n = 0.105;
    const M = n * t;

    let E = M;
    for (let i = 0; i < 7; i++) {
      E -= (E - eccentricity * Math.sin(E) - M) / (1.0 - eccentricity * Math.cos(E));
    }

    const cosE = Math.cos(E);
    const sinE = Math.sin(E);
    const trueAnomaly = 2.0 * Math.atan2(Math.sqrt(1 + eccentricity) * Math.sin(E * 0.5), Math.sqrt(1 - eccentricity) * Math.cos(E * 0.5));
    const precession = 0.22 * t;
    const orbitalAngle = trueAnomaly + precession;
    const radius = p / (1.0 + eccentricity * Math.cos(trueAnomaly));

    const x = radius * Math.cos(orbitalAngle);
    const z = radius * Math.sin(orbitalAngle);

    const dThetaDt = n * Math.sqrt(1.0 - eccentricity * eccentricity) / Math.max(1e-3, (1.0 - eccentricity * cosE));
    const vxOrb = -radius * Math.sin(orbitalAngle) * dThetaDt;
    const vzOrb = radius * Math.cos(orbitalAngle) * dThetaDt;
    const tangent = normalize([vxOrb, 0.0, vzOrb]);

    const skimHeight = 0.25 * Math.sin(orbitalAngle * 0.5);
    const camPos = [x, skimHeight, z];
    const vy = 0.25 * 0.5 * Math.cos(orbitalAngle * 0.5) * dThetaDt;
    const velLen = Math.hypot(vxOrb, vy, vzOrb) || 1e-6;
    const betaKepler = Math.min(1.0 / Math.sqrt(Math.max(radius, 3.0)), 0.85);
    const camVel = [
      vxOrb / velLen * betaKepler,
      vy / velLen * betaKepler,
      vzOrb / velLen * betaKepler,
    ];
    const forward = tangent;

    const worldUp = [0, 1, 0];
    let right = cross(forward, worldUp);
    if (Math.hypot(right[0], right[1], right[2]) < 1e-4) right = [1, 0, 0];
    right = normalize(right);
    const up = normalize(cross(right, forward));

    return { camPos, right, up, forward, camVel };
  }

  function getManualCameraBasis() {
    const yaw = camYaw * Math.PI / 180;
    const pitch = camPitch * Math.PI / 180;
    const cp = Math.cos(pitch);
    const sp = Math.sin(pitch);
    const cy = Math.cos(yaw);
    const sy = Math.sin(yaw);

    const camPos = [
      camDist * cp * sy,
      camDist * sp,
      camDist * cp * cy,
    ];

    const forward = normalize([-camPos[0], -camPos[1], -camPos[2]]);
    const worldUp = [0, 1, 0];
    let right = cross(forward, worldUp);
    if (Math.hypot(right[0], right[1], right[2]) < 1e-4) {
      right = [1, 0, 0];
    }
    right = normalize(right);
    const up = normalize(cross(right, forward));

    return { camPos, right, up, forward, camVel: [0, 0, 0] };
  }

  incSlider.addEventListener('input', (e) => {
    incl = parseFloat(e.target.value) * Math.PI / 180;
    updateReadouts();
  });

  routSlider.addEventListener('input', (e) => {
    rout = parseFloat(e.target.value);
    updateReadouts();
  });

  rotSlider.addEventListener('input', (e) => {
    speed = parseFloat(e.target.value);
    updateReadouts();
  });

  yawSlider.addEventListener('input', (e) => {
    camYaw = parseFloat(e.target.value);
    updateReadouts();
  });

  pitchSlider.addEventListener('input', (e) => {
    camPitch = parseFloat(e.target.value);
    updateReadouts();
  });

  distSlider.addEventListener('input', (e) => {
    camDist = parseFloat(e.target.value);
    updateReadouts();
  });

  thickSlider.addEventListener('input', (e) => {
    diskHalfThickness = parseFloat(e.target.value);
    updateReadouts();
  });

  bendSlider.addEventListener('input', (e) => {
    bend = parseFloat(e.target.value);
    updateReadouts();
  });

  spinSlider.addEventListener('input', (e) => {
    spin = parseFloat(e.target.value);
    updateReadouts();
  });

  progradeToggle.addEventListener('change', (e) => {
    diskPrograde = e.target.checked;
    updateReadouts();
  });

  qualitySlider.addEventListener('input', (e) => {
    quality = parseFloat(e.target.value);
    updateReadouts();
  });

  autoLodToggle.addEventListener('change', (e) => {
    autoLodEnabled = e.target.checked;
    updateReadouts();
  });

  lodMinSlider.addEventListener('input', (e) => {
    lodMin = parseFloat(e.target.value);
    if (lodMin > lodMax - 0.01) lodMax = clamp(lodMin + 0.01, 0.25, 1.0);
    updateReadouts();
  });

  lodMaxSlider.addEventListener('input', (e) => {
    lodMax = parseFloat(e.target.value);
    if (lodMax < lodMin + 0.01) lodMin = clamp(lodMax - 0.01, 0.25, 1.0);
    updateReadouts();
  });

  maxDiskCrossingsSlider.addEventListener('input', (e) => {
    maxDiskCrossings = parseInt(e.target.value, 10);
    updateReadouts();
  });

  followOrbitToggle.addEventListener('change', (e) => {
    followOrbit = e.target.checked;
    yawSlider.disabled = followOrbit;
    pitchSlider.disabled = followOrbit;
    distSlider.disabled = followOrbit;
    setModeBadge(followOrbit
      ? "▸ CAMÉRA ORBITALE — géodésique excentrique (style Interstellar)"
      : "▸ CAMÉRA MANUELLE — Schwarzschild/Kerr-lite",
      followOrbit ? "geodesic" : "");
  });

  emissionModeToggle.addEventListener('change', (e) => {
    ntEmissionMode = e.target.checked;
    updateReadouts();
  });

  modeNormal.addEventListener('click', () => setBloomMode(false));
  modeBloom.addEventListener('click', () => setBloomMode(true));

  document.addEventListener('keydown', (e) => {
    const key = e.key.toLowerCase();
    if (key === 'a') setModeBadge('▸ ANALYTIQUE — Luminet Eq.A3');
    if (key === 'g') setModeBadge('▸ GÉODÉSIQUE RK4 — Schwarzschild/Kerr-lite', 'geodesic');
    if (key === 'l') setModeBadge('▸ LUT — mode expérimental', 'lut');
  });

  let dragging = false;
  function pointerToParams(clientX, clientY) {
    const xNorm = clientX / window.innerWidth;
    const yNorm = clientY / window.innerHeight;
    incl = (5 + yNorm * 85) * Math.PI / 180;
    rout = 10 + xNorm * 50;
    camYaw = -180 + xNorm * 360;
    camPitch = clamp(70 - yNorm * 140, -75, 75);
    updateReadouts();
  }

  window.addEventListener('pointerdown', (e) => {
    dragging = true;
    crosshair.style.display = 'block';
    pointerToParams(e.clientX, e.clientY);
  });

  window.addEventListener('pointermove', (e) => {
    crosshair.style.left = `${e.clientX}px`;
    crosshair.style.top = `${e.clientY}px`;
    if (dragging) pointerToParams(e.clientX, e.clientY);
  });

  window.addEventListener('pointerup', () => {
    dragging = false;
    crosshair.style.display = 'none';
  });

  function resize() {
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.floor(canvas.clientWidth * dpr);
    canvas.height = Math.floor(canvas.clientHeight * dpr);
    gl.viewport(0, 0, canvas.width, canvas.height);
  }
  window.addEventListener('resize', resize);
  resize();
  updateReadouts();
  setBloomMode(false);
  followOrbitToggle.checked = true;
  followOrbitToggle.dispatchEvent(new Event('change'));

  let frames = 0;
  let fpsTimer = 0;
  const fpsSpan = document.getElementById('fps');

  function render(now) {
    now *= 0.001;
    frames++;
    if (now - fpsTimer >= 1.0) {
      fpsSpan.textContent = `${Math.round(frames / (now - fpsTimer))} fps`;
      frames = 0;
      fpsTimer = now;
    }

    const orbitTime = now * speed;
    const basis = followOrbit ? getSkimmingOrbitBasis(orbitTime) : getManualCameraBasis();

    const rawCamSpeed = Math.hypot(basis.camVel[0], basis.camVel[1], basis.camVel[2]);
    const schmittBand = 0.02;
    let targetSpeed = filteredCamSpeed;
    if (rawCamSpeed > filteredCamSpeed + schmittBand) {
      targetSpeed = rawCamSpeed - schmittBand;
    } else if (rawCamSpeed < filteredCamSpeed - schmittBand) {
      targetSpeed = rawCamSpeed + schmittBand;
    }
    const dt = Math.min(Math.max(now - (render.lastNow || now), 0.0), 0.2);
    const riseTau = 0.22;
    const fallTau = 0.38;
    const tau = targetSpeed >= filteredCamSpeed ? riseTau : fallTau;
    const alpha = 1.0 - Math.exp(-dt / Math.max(tau, 1e-4));
    filteredCamSpeed += (targetSpeed - filteredCamSpeed) * alpha;
    render.lastNow = now;

    gl.uniform2f(uRes, canvas.width, canvas.height);
    gl.uniform1f(uTime, now);
    gl.uniform1f(uIncl, incl);
    gl.uniform1f(uRout, rout);
    gl.uniform1f(uSpeed, speed);
    gl.uniform1f(uBloom, bloomMode);
    gl.uniform3f(uCamPos, basis.camPos[0], basis.camPos[1], basis.camPos[2]);
    gl.uniform3f(uCamRight, basis.right[0], basis.right[1], basis.right[2]);
    gl.uniform3f(uCamUp, basis.up[0], basis.up[1], basis.up[2]);
    gl.uniform3f(uCamForward, basis.forward[0], basis.forward[1], basis.forward[2]);
    gl.uniform1f(uDiskHalfThickness, diskHalfThickness);
    gl.uniform1f(uBend, bend);
    gl.uniform1f(uSpin, spin);
    gl.uniform1f(uQuality, quality);
    gl.uniform1f(uFovX, fovX);
    gl.uniform3f(uCamVel, basis.camVel[0], basis.camVel[1], basis.camVel[2]);
    gl.uniform1f(uFogStrength, fogStrength);
    gl.uniform1f(uDiskSense, diskPrograde ? 1.0 : 0.0);
    gl.uniform1f(uEmissionMode, ntEmissionMode ? 1.0 : 0.0);
    gl.uniform1f(uMaxDiskCrossings, maxDiskCrossings);
    gl.uniform1f(uLodAuto, autoLodEnabled ? 1.0 : 0.0);
    gl.uniform1f(uLodMin, lodMin);
    gl.uniform1f(uLodMax, lodMax);
    gl.uniform1f(uCamSpeedFiltered, filteredCamSpeed);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    requestAnimationFrame(render);
  }
  requestAnimationFrame(render);
})();
