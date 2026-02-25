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

  #define PI 3.141592653589793
  #define ISCO 6.0
  #define BC 5.196152
  #define HORIZON 2.0
  #define MAX_STEPS 300
  #define SKY_RADIUS 180.0

  float hash(vec2 p){ return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453); }
  float redshift(float r,float phi,float incl){ float v=1.0/sqrt(max(r-3.0,0.001)); return max((1.0+v*cos(phi)*sin(incl))/sqrt(max(1.0-3.0/r,0.001)),0.001); }
  float emissivity(float r,float phi,float t){
    if(r<ISCO || r>iRout) return 0.0;
    float base=(1.0-sqrt(ISCO/r))*pow(r,-3.0);
    float spin=fract(phi/(2.0*PI)+t*0.1*iSpeed);
    float flare=pow(sin(spin*20.0+r*0.5)*0.5+0.5,4.0)*0.3;
    return base*(1.0+flare);
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
  vec3 aces(vec3 x){ const float a=2.51,b=0.03,c=2.43,d=0.59,e=0.14; return clamp((x*(a*x+b))/(x*(c*x+d)+e),0.0,1.0); }
  vec3 bloom(vec3 col,float intensity){ float l=dot(col,vec3(0.2126,0.7152,0.0722)); if(l>0.8) col+=vec3(0.4,0.3,0.2)*(l-0.8)*2.0*intensity; return col; }

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

  vec3 skyColor(vec3 dir){
    vec3 d = normalize(dir);
    float phi = atan(d.z, d.x);
    float theta = acos(clamp(d.y, -1.0, 1.0));

    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);
    vec2 g = floor(uv * vec2(180.0, 90.0));

    float cluster = hash(g);
    float stars = step(0.994, cluster) * (0.25 + 0.75 * hash(g * 0.73 + 9.1));

    float milkyBand = exp(-pow((uv.y - 0.5) * 7.0, 2.0));
    float dustBand = hash(g * vec2(1.3, 2.1) + 17.0) * milkyBand;

    vec3 nebula = mix(vec3(0.03, 0.035, 0.055), vec3(0.08, 0.06, 0.04), dustBand * 0.8);
    vec3 starTint = mix(vec3(0.85, 0.9, 1.0), vec3(1.0, 0.88, 0.72), hash(g + 3.7));

    return nebula + starTint * stars * (0.7 + 0.3 * milkyBand);
  }

  bool traceScene(vec3 ro, vec3 rd, out vec3 hit, out vec3 hitDir, out vec3 skyDir, out float bMin, out bool swallowed){
    vec3 pos = ro;
    vec3 dir = rd;
    float L2 = dot(cross(ro, rd), cross(ro, rd));
    float prevY = pos.y;
    bMin = 1e9;
    swallowed = false;
    skyDir = rd;

    for (int i = 0; i < MAX_STEPS; i++) {
      if (float(i) > mix(90.0, float(MAX_STEPS), iQuality)) break;

      float r = length(pos);
      bMin = min(bMin, length(cross(pos, dir)));
      if (r < HORIZON) {
        swallowed = true;
        return false;
      }
      if (r > SKY_RADIUS) {
        skyDir = normalize(pos);
        return false;
      }

      float stepFar = mix(0.45, 1.4, clamp((r - 8.0) / 55.0, 0.0, 1.0));
      float h = stepFar / mix(1.5, 0.7, iQuality);
      h *= mix(1.0, 0.72, clamp(10.0 / max(r, 0.1), 0.0, 1.0));

      vec3 prevPos = pos;
      vec3 prevDir = dir;
      rk4Step(pos, dir, h, L2);

      if (sign(prevY) != sign(pos.y)) {
        float t = prevY / (prevY - pos.y);
        vec3 crossPos = mix(prevPos, pos, t);
        if (abs(crossPos.y) <= iDiskHalfThickness && length(crossPos.xz) > ISCO && length(crossPos.xz) < iRout) {
          hit = crossPos;
          hitDir = normalize(mix(prevDir, dir, t));
          return true;
        }
      }
      prevY = pos.y;
    }
    skyDir = normalize(pos);
    return false;
  }

  void main(){
    vec2 p=(gl_FragCoord.xy-0.5*iResolution.xy)/iResolution.y;
    float focal = 0.5 * (iResolution.x / iResolution.y) / tan(radians(iFovX) * 0.5);
    vec3 rdCam = normalize(vec3(p.x, p.y, focal));
    vec3 rayDir = normalize(iCamRight * rdCam.x + iCamUp * rdCam.y + iCamForward * rdCam.z);

    vec3 hit = vec3(0.0);
    vec3 hitDir = vec3(0.0);
    vec3 skyDir = vec3(0.0);
    float bMin;
    bool swallowed = false;
    bool seenDisk = traceScene(iCamPos, rayDir, hit, hitDir, skyDir, bMin, swallowed);

    vec3 bg = skyColor(skyDir);

    if (swallowed) {
      float photonRing = exp(-pow((bMin - BC) * 3.5, 2.0));
      bg = vec3(0.0);
      bg += vec3(0.95, 0.66, 0.35) * photonRing * 0.45;
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

    if(seenDisk){
      float rHit = length(hit.xz);
      float alpha = atan(hit.z, hit.x);

      float z1=redshift(rHit,alpha,iIncl);
      float flux1=emissivity(rHit,alpha,iTime)*pow(1.0/z1,4.0);
      flux1*=dust(vec2(rHit,alpha),iTime)*170.0;
      vec3 c1=tempColor(clamp(0.5+(1.0/z1-1.0)*1.2,0.0,1.0))*flux1;

      float rGhost = rHit + 10.0 * exp(-abs(bMin - BC) * 1.7);
      float z2=redshift(rGhost,alpha+PI,iIncl);
      float flux2=emissivity(rGhost,alpha+PI,iTime)*pow(1.0/z2,4.0);
      flux2*=dust(vec2(rGhost,alpha+PI),iTime)*40.0;
      vec3 c2=tempColor(clamp(0.5+(1.0/z2-1.0)*1.2,0.0,1.0))*flux2*0.28;

      float thicknessFade = smoothstep(iDiskHalfThickness, 0.0, abs(hit.y));
      float dopplerAniso = 0.75 + 0.25 * clamp(dot(normalize(vec3(-hit.z,0.0,hit.x)), -hitDir), -1.0, 1.0);
      col += (c1 + c2) * thicknessFade * dopplerAniso;
    }

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
  const qualitySlider = document.getElementById('qualitySlider');
  const qualitySliderVal = document.getElementById('qualitySliderVal');
  const followOrbitToggle = document.getElementById('followOrbitToggle');
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
  let quality = 0.72;
  const fovX = 90.0;

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
    qualitySlider.value = quality.toFixed(2);
    qualitySliderVal.textContent = quality.toFixed(2);
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
    const forward = tangent;

    const worldUp = [0, 1, 0];
    let right = cross(forward, worldUp);
    if (Math.hypot(right[0], right[1], right[2]) < 1e-4) right = [1, 0, 0];
    right = normalize(right);
    const up = normalize(cross(right, forward));

    return { camPos, right, up, forward };
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

    return { camPos, right, up, forward };
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

  qualitySlider.addEventListener('input', (e) => {
    quality = parseFloat(e.target.value);
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
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    requestAnimationFrame(render);
  }
  requestAnimationFrame(render);
})();
