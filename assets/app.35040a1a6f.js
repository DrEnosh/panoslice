(() => {
'use strict';

/* =====================================================================
   Shared reconstruction code. It runs on the main thread and, through
   Function#toString, inside a Web Worker, so it must stay self-contained.
   World frame (mm): +x = patient left, +y = anterior, +z = superior.
   ===================================================================== */
function shared() {
  const clamp = (v, a, b) => v < a ? a : v > b ? b : v;
  const lerp = (a, b, t) => a + (b - a) * t;
  const smoothstep = (a, b, x) => { const t = clamp((x - a) / (b - a), 0, 1); return t * t * (3 - 2 * t); };

  // Half arch, midline incisors -> condyle (patient left). Mirrored for the right.
  const ARCH_PTS = [[0, 0], [9, -2], [16, -7], [21, -14], [24, -21], [26.5, -30], [28, -40], [29, -49], [32, -60], [36, -72]];
  // Mean adult tooth dimensions (mm): mesiodistal width, crown height, root length, visible roots.
  const UPPER = [{w:8.5,ch:10.5,rl:13,r:1},{w:6.5,ch:9,rl:13,r:1},{w:7.5,ch:10,rl:17,r:1},{w:7,ch:8.5,rl:14,r:2},{w:6.5,ch:8,rl:14,r:1},{w:10,ch:7.5,rl:12.5,r:3},{w:9,ch:7,rl:12,r:3},{w:8.5,ch:6.5,rl:11,r:2}];
  const LOWER = [{w:5.5,ch:9,rl:12.5,r:1},{w:6,ch:9.5,rl:14,r:1},{w:7,ch:11,rl:16,r:1},{w:7,ch:8.5,rl:14,r:1},{w:7.5,ch:8,rl:14.5,r:1},{w:11,ch:7.5,rl:14,r:2},{w:10.5,ch:7,rl:13,r:2},{w:10,ch:7,rl:11,r:2}];
  // Bucco-lingual thickness (mm) against fraction of the half-arch length.
  const THICK = [[0, 7.5], [0.14, 8], [0.23, 9], [0.35, 9.5], [0.48, 11], [0.72, 11], [0.8, 8.5], [1, 7]];
  function thickAt(f) {
    f = Math.abs(f);
    for (let i = 1; i < THICK.length; i++) if (f <= THICK[i][0]) {
      const [a, va] = THICK[i - 1], [b, vb] = THICK[i];
      return lerp(va, vb, (f - a) / (b - a));
    }
    return THICK[THICK.length - 1][1];
  }
  function buildArch(wS, dS) {
    const half = ARCH_PTS.map(([x, y]) => [x * wS, y * dS]);
    const pts = [];
    for (let i = half.length - 1; i > 0; i--) pts.push([-half[i][0], half[i][1]]);
    for (const p of half) pts.push(p);
    const dense = [], m = pts.length;
    for (let i = 0; i < m - 1; i++) {
      const p0 = pts[Math.max(0, i - 1)], p1 = pts[i], p2 = pts[i + 1], p3 = pts[Math.min(m - 1, i + 2)];
      for (let k = 0; k < 32; k++) {
        const t = k / 32, t2 = t * t, t3 = t2 * t;
        const f = (a, b, c, d) => 0.5 * (2 * b + (-a + c) * t + (2 * a - 5 * b + 4 * c - d) * t2 + (-a + 3 * b - 3 * c + d) * t3);
        dense.push([f(p0[0], p1[0], p2[0], p3[0]), f(p0[1], p1[1], p2[1], p3[1])]);
      }
    }
    dense.push(pts[m - 1]);
    const cum = [0];
    for (let i = 1; i < dense.length; i++) cum.push(cum[i - 1] + Math.hypot(dense[i][0] - dense[i - 1][0], dense[i][1] - dense[i - 1][1]));
    const total = cum[cum.length - 1], N = 1024;
    const x = new Float32Array(N), y = new Float32Array(N), s = new Float32Array(N);
    const tx = new Float32Array(N), ty = new Float32Array(N), nx = new Float32Array(N), ny = new Float32Array(N);
    let j = 0;
    for (let i = 0; i < N; i++) {
      const target = total * i / (N - 1);
      while (j < cum.length - 2 && cum[j + 1] < target) j++;
      const u = (target - cum[j]) / ((cum[j + 1] - cum[j]) || 1);
      x[i] = lerp(dense[j][0], dense[j + 1][0], u);
      y[i] = lerp(dense[j][1], dense[j + 1][1], u);
      s[i] = target - total / 2;
    }
    for (let i = 0; i < N; i++) {
      const a = Math.max(0, i - 1), b = Math.min(N - 1, i + 1);
      let dx = x[b] - x[a], dy = y[b] - y[a]; const l = Math.hypot(dx, dy) || 1; dx /= l; dy /= l;
      tx[i] = dx; ty[i] = dy; nx[i] = -dy; ny[i] = dx; // outward (buccal) normal
    }
    return { N, x, y, s, tx, ty, nx, ny, L: total / 2 };
  }
  const REF_L = buildArch(1, 1).L;
  for (const T of [UPPER, LOWER]) { let acc = 0; for (const t of T) { t.start = acc / REF_L; t.c = acc + t.w / 2; acc += t.w; t.f = t.c / REF_L; t.fEnd = acc / REF_L; } }

  function nearestIdx(A, x, y) {
    const X = A.x, Y = A.y, N = A.N; let best = Infinity, bk = 0;
    for (let k = 0; k < N; k += 8) { const dx = x - X[k], dy = y - Y[k], d = dx * dx + dy * dy; if (d < best) { best = d; bk = k; } }
    const k1 = Math.min(N - 1, bk + 8);
    for (let k = Math.max(0, bk - 8); k <= k1; k++) { const dx = x - X[k], dy = y - Y[k], d = dx * dx + dy * dy; if (d < best) { best = d; bk = k; } }
    return bk;
  }
  const idxAtS = (A, s) => clamp(Math.round((s + A.L) / (2 * A.L) * (A.N - 1)), 0, A.N - 1);
  function tableAt(T, af, key) {
    if (af <= T[0].f) return T[0][key];
    for (let i = 1; i < T.length; i++) if (af <= T[i].f) return lerp(T[i - 1][key], T[i][key], (af - T[i - 1].f) / (T[i].f - T[i - 1].f));
    return T[T.length - 1][key];
  }
  function toothIndexAt(T, af) { for (let i = 0; i < T.length; i++) if (af < T[i].fEnd) return i; return T.length - 1; }
  // The occlusal plane on the OPG is a "smile": highest at the ends.
  const occFracAt = (P, u) => { const n = (u - P.mid) / (P.span / 2); return P.occl - P.curve * n * n; };

  function boxPass(src, dst, w, h, r, horiz) {
    const inv = 1 / (2 * r + 1);
    if (horiz) {
      for (let y = 0; y < h; y++) {
        const o = y * w; let acc = 0;
        for (let k = -r; k <= r; k++) acc += src[o + (k < 0 ? 0 : k >= w ? w - 1 : k)];
        for (let x = 0; x < w; x++) {
          dst[o + x] = acc * inv;
          const a = x + r + 1, b = x - r;
          acc += src[o + (a >= w ? w - 1 : a)] - src[o + (b < 0 ? 0 : b)];
        }
      }
    } else {
      for (let x = 0; x < w; x++) {
        let acc = 0;
        for (let k = -r; k <= r; k++) acc += src[(k < 0 ? 0 : k >= h ? h - 1 : k) * w + x];
        for (let y = 0; y < h; y++) {
          dst[y * w + x] = acc * inv;
          const a = y + r + 1, b = y - r;
          acc += src[(a >= h ? h - 1 : a) * w + x] - src[(b < 0 ? 0 : b) * w + x];
        }
      }
    }
  }
  const radiusFor = s => Math.max(1, Math.round((Math.sqrt(1 + 4 * s * s) - 1) / 2));
  function blurXY(src, w, h, sx, sy) {
    const out = Float32Array.from(src), tmp = new Float32Array(src.length);
    const rx = sx > 0.35 ? radiusFor(sx) : 0, ry = sy > 0.35 ? radiusFor(sy) : 0;
    for (let p = 0; p < 3; p++) {
      if (rx) { boxPass(out, tmp, w, h, rx, true); out.set(tmp); }
      if (ry) { boxPass(out, tmp, w, h, ry, false); out.set(tmp); }
    }
    return out;
  }
  const gaussBlur = (src, w, h, s) => blurXY(src, w, h, s, s);
  function percentileNormalize(g) {
    const n = g.length, hist = new Uint32Array(1024);
    for (let i = 0; i < n; i++) hist[Math.min(1023, Math.max(0, (g[i] * 1023) | 0))]++;
    const pick = q => { let acc = 0; const t = q * n; for (let i = 0; i < 1024; i++) { acc += hist[i]; if (acc >= t) return i / 1023; } return 1; };
    const lo = pick(0.005), hi = Math.max(lo + 0.02, pick(0.997));
    const out = new Float32Array(n), sc = 1 / (hi - lo);
    for (let i = 0; i < n; i++) out[i] = clamp((g[i] - lo) * sc, 0, 1);
    return out;
  }
  function rangeH(src, w, h, r, isMax) {
    const out = new Float32Array(src.length);
    for (let y = 0; y < h; y++) {
      const o = y * w;
      for (let x = 0; x < w; x++) {
        let m = src[o + x];
        const a = Math.max(0, x - r), b = Math.min(w - 1, x + r);
        if (isMax) { for (let k = a; k <= b; k++) if (src[o + k] > m) m = src[o + k]; }
        else { for (let k = a; k <= b; k++) if (src[o + k] < m) m = src[o + k]; }
        out[o + x] = m;
      }
    }
    return out;
  }
  // 2D "toothness": teeth are brighter than the bone or air right beside them.
  function toothMap(norm, W, H, ppm) {
    const Ib = gaussBlur(norm, W, H, 0.8);
    const bg = blurXY(norm, W, H, 4 * ppm, 1.5 * ppm);
    const c = new Float32Array(W * H);
    for (let i = 0; i < c.length; i++) c[i] = Math.max(smoothstep(0.015, 0.09, Ib[i] - bg[i]), smoothstep(0.66, 0.84, Ib[i]));
    const r = Math.max(1, Math.round(0.9 * ppm));          // close pulp canals and PDL slits
    return gaussBlur(rangeH(rangeH(c, W, H, r, true), W, H, r, false), W, H, 0.8);
  }
  function blur3(a, nx, ny, nz) {
    const line = new Float32Array(Math.max(nx, ny, nz));
    const dims = [nx, ny, nz], strides = [4, nx * 4, nx * ny * 4];
    for (let ax = 0; ax < 3; ax++) {
      const n = dims[ax], st = strides[ax];
      const o1 = ax === 0 ? ny : nx, o2 = ax === 2 ? ny : nz;
      for (let b = 0; b < o2; b++) for (let c = 0; c < o1; c++) {
        const start = ax === 0 ? (c * nx + b * nx * ny) * 4 : ax === 1 ? (c + b * nx * ny) * 4 : (c + b * nx) * 4;
        for (let ch = 0; ch < 4; ch++) {
          let any = 0;
          for (let i = 0; i < n; i++) { const v = a[start + ch + i * st]; line[i] = v; any |= v; }
          if (!any) continue;
          for (let i = 0; i < n; i++) {
            const p = i > 0 ? line[i - 1] : line[i], q = i < n - 1 ? line[i + 1] : line[i];
            a[start + ch + i * st] = (p + 2 * line[i] + q + 2) >> 2;
          }
        }
      }
    }
  }
  function bil(M, a, b, c, d, fu, fv) { return (M[a] * (1 - fu) + M[b] * fu) * (1 - fv) + (M[c] * (1 - fu) + M[d] * fu) * fv; }

  /* Builds an RGBA volume: R = hard-tissue density, G = toothness, B = gum, A = metal. */
  function buildVolume(norm, raw, W, H, P, cache, teeth) {
    const t0 = performance.now();
    const A = buildArch(P.archW, P.archD), L = A.L, N = A.N;
    const mmPerPx = (2 * L) / (P.span * W), ppm = 1 / mmPerPx;
    const Hmm = H * mmPerPx * P.vscale;
    let maxX = 0, minY = 0, maxY = -1e9;
    for (let i = 0; i < N; i++) { maxX = Math.max(maxX, Math.abs(A.x[i])); minY = Math.min(minY, A.y[i]); maxY = Math.max(maxY, A.y[i]); }
    const lat = 11 * P.thick / 2 + 6;
    const bmin = [-maxX - lat, minY - 6, 0];
    const size0 = [2 * (maxX + lat), maxY + 16 - bmin[1], Hmm];
    const target = P.res === 'high' ? 240 : P.res === 'lite' ? 128 : 176;
    const vox = Math.max(size0[0], size0[1], size0[2]) / target;
    const nx = Math.ceil(size0[0] / vox), ny = Math.ceil(size0[1] / vox), nz = Math.ceil(size0[2] / vox);
    const bmax = [bmin[0] + nx * vox, bmin[1] + ny * vox, bmin[2] + nz * vox];

    const sI = Math.hypot(P.smooth, 0.45 * vox * ppm);
    const sS = sI + 1.5 + 0.5 * ppm, sG = 2.5 * ppm;
    const key = [W, H, sI.toFixed(2), sS.toFixed(2), sG.toFixed(2), ppm.toFixed(3)].join('|');
    if (cache.key !== key || cache.src !== norm) {
      cache.I = gaussBlur(norm, W, H, sI); cache.Is = gaussBlur(norm, W, H, sS); cache.Ig = blurXY(norm, W, H, sG * 2, sG * 0.5);
      cache.Rw = gaussBlur(raw, W, H, 0.8); cache.key = key; cache.src = norm; cache.tmKey = null;
    }
    const tmKey = teeth ? 'model' : 'contrast|' + ppm.toFixed(3);
    if (cache.tmKey !== tmKey || cache.tmFor !== (teeth ? teeth.prob : norm)) {
      cache.Tm = teeth ? teethToMap(teeth, W, H) : toothMap(norm, W, H, ppm);
      cache.tmKey = tmKey; cache.tmFor = teeth ? teeth.prob : norm;
    }
    const I = cache.I, Is = cache.Is, Ig = cache.Ig, Tm = cache.Tm, Rw = cache.Rw;

    // Per-arch-sample lookups
    const uPx = new Float32Array(N), zOcc = new Float32Array(N), chU = new Float32Array(N), chL = new Float32Array(N);
    const lenU = new Float32Array(N), lenL = new Float32Array(N), mU = new Float32Array(N), mL = new Float32Array(N);
    const wa = new Float32Array(N), TK = new Float32Array(N), bas = new Float32Array(N);
    const dentU = new Float32Array(N), dentL = new Float32Array(N), gumK = new Float32Array(N), capU = new Float32Array(N);
    const fEU = UPPER[7].fEnd, fEL = LOWER[7].fEnd, fEM = Math.max(fEU, fEL);
    for (let k = 0; k < N; k++) {
      const s = A.s[k], af = Math.abs(s) / L, u = P.mid + P.span * s / (2 * L);
      uPx[k] = u * W - 0.5;
      zOcc[k] = (1 - occFracAt(P, u)) * Hmm;
      chU[k] = tableAt(UPPER, af, 'ch'); chL[k] = tableAt(LOWER, af, 'ch');
      lenU[k] = chU[k] + tableAt(UPPER, af, 'rl') + 1.5; lenL[k] = chL[k] + tableAt(LOWER, af, 'rl') + 1.5;
      const iu = toothIndexAt(UPPER, af), il = toothIndexAt(LOWER, af);
      const phU = clamp((af - UPPER[iu].start) / (UPPER[iu].fEnd - UPPER[iu].start), 0, 1);
      const phL = clamp((af - LOWER[il].start) / (LOWER[il].fEnd - LOWER[il].start), 0, 1);
      const scU = af > fEU ? 0 : 2.4 * Math.pow(1 - Math.sin(Math.PI * phU), 1.6);
      const scL = af > fEL ? 0 : 2.2 * Math.pow(1 - Math.sin(Math.PI * phL), 1.6);
      mU[k] = chU[k] - 1.1 - scU; mL[k] = chL[k] - 1.1 - scL;
      wa[k] = 1 - smoothstep(0.1, 0.28, af);
      TK[k] = thickAt(af) * P.thick; bas[k] = lerp(2, 6, wa[k]) * P.thick;
      dentU[k] = 1 - smoothstep(fEU + 0.015, fEU + 0.05, af); dentL[k] = 1 - smoothstep(fEL + 0.015, fEL + 0.05, af);
      gumK[k] = 1 - smoothstep(fEM + 0.05, fEM + 0.11, af);
      capU[k] = 1 - smoothstep(fEU + 0.02, fEU + 0.09, af);   // round off the maxilla above the apices, but not the ramus
    }

    const nvox = nx * ny * nz, rgba = new Uint8Array(nvox * 4);
    const cmin = [Infinity, Infinity, Infinity], cmax = [-Infinity, -Infinity, -Infinity];
    const edge = Math.max(vox * 1.1, 0.5), fl = P.floor, sk = P.shell, gam = P.round, stride4 = nx * ny * 4;
    const rowF = new Float32Array(nz), zc = new Float32Array(nz);
    for (let k = 0; k < nz; k++) { zc[k] = bmin[2] + (k + 0.5) * vox; rowF[k] = (1 - zc[k] / Hmm) * H - 0.5; }
    for (let j = 0; j < ny; j++) {
      const y = bmin[1] + (j + 0.5) * vox;
      for (let i = 0; i < nx; i++) {
        const x = bmin[0] + (i + 0.5) * vox;
        const k = nearestIdx(A, x, y);
        if (k <= 1 || k >= N - 2) continue;
        const d = (x - A.x[k]) * A.nx[k] + (y - A.y[k]) * A.ny[k];
        const wak = wa[k];
        if (Math.abs(d) > TK[k] / 2 + bas[k] / 2 + 4.5 + 10 * wak) continue;
        const u = uPx[k]; if (u < 0 || u > W - 1) continue;
        const u0 = u | 0, fu = u - u0, u1 = Math.min(u0 + 1, W - 1);
        const base = (i + nx * j) * 4, zo = zOcc[k];
        for (let kz = 0; kz < nz; kz++) {
          const v = rowF[kz]; if (v < 0 || v > H - 1) continue;
          const dz = zc[kz] - zo;
          const up = dz <= -1.5 ? 0 : dz >= 1.5 ? 1 : smoothstep(-1.5, 1.5, dz);
          // Upper arch sits a little buccal to the lower; incisors procline; the chin comes forward.
          let off = -0.6 + 1.5 * up;
          if (wak > 0) off += wak * (up * (-0.33 * clamp(dz, 0, 24)) + (1 - up) * (-0.36 * clamp(-dz, 0, 22) + 0.45 * clamp(-dz - 24, 0, 14)));
          const ad = Math.abs(d - off);
          let hwMax = TK[k] / 2;
          if (dz < -22) hwMax += bas[k] / 2 * smoothstep(22, 32, -dz);
          else if (dz > 26) hwMax *= 1 - capU[k] * smoothstep(26, 34, dz);
          if (ad > hwMax + 3) continue;
          const v0 = v | 0, fv = v - v0, r0 = v0 * W, r1 = Math.min(v0 + 1, H - 1) * W;
          const a00 = r0 + u0, a01 = r0 + u1, a10 = r1 + u0, a11 = r1 + u1;
          const q = (bil(Is, a00, a01, a10, a11, fu, fv) - fl) / (1 - fl);
          let hw = q > 0 ? hwMax * (0.22 + 0.78 * Math.pow(q > 1 ? 1 : q, gam)) * smoothstep(0, 0.06, q) : 0;
          let R = 0, G = 0, B = 0, M = 0, fill = 0;
          {
            const adz = dz < 0 ? -dz : dz, ch = dz >= 0 ? chU[k] : chL[k], len = dz >= 0 ? lenU[k] : lenL[k];
            const alv = smoothstep(ch - 0.5, ch + 1.5, adz) * (dz >= 0 ? (1 - smoothstep(24, 30, adz)) * capU[k] : gumK[k]);
            if (alv > 0) {
              const qa = (bil(Ig, a00, a01, a10, a11, fu, fv) - fl * 0.75) / (1 - fl * 0.75);
              if (qa > 0) { const hwA = hwMax * 0.72 * alv * Math.min(1, qa * 4); fill = 1; hw = hw > hwA ? lerp(hw, hwA, 0.65 * alv) : hwA; }
            }
          }
          if (hw > 0) {
            let w = (hw - ad) / edge + 0.5;
            if (w > 0) {
              if (w > 1) w = 1;
              const val = bil(I, a00, a01, a10, a11, fu, fv);
              let sh = (ad - 0.45 * hw) / (0.55 * hw); sh = sh < 0 ? 0 : sh > 1 ? 1 : sh; sh = sh * sh * (3 - 2 * sh);
              R = val * ((1 - 0.45 * sk) + 0.9 * sk * sh) * w;
              if (fill) R = Math.max(R, (0.48 + 0.1 * sh) * w);
              const rv = bil(Rw, a00, a01, a10, a11, fu, fv);
              M = smoothstep(0.93, 0.97, rv) * w;
              const tm = bil(Tm, a00, a01, a10, a11, fu, fv);
              if (tm > 0.02) {
                const adz = dz < 0 ? -dz : dz;
                const len = dz >= 0 ? lenU[k] : lenL[k], ch = dz >= 0 ? chU[k] : chL[k];
                const band = (1 - smoothstep(len - 2, len + 2, adz)) * (dz >= 0 ? dentU[k] : dentL[k]);
                const kt = 0.8 + 0.2 * (1 - smoothstep(ch - 1, ch + 1, adz));
                let wt = (hw * kt - ad) / edge + 0.5; wt = wt < 0 ? 0 : wt > 1 ? 1 : wt;
                G = Math.max(tm * band * wt, M);
              } else G = M;
            }
          }
          const zg = dz >= 0
            ? smoothstep(mU[k] - 0.4, mU[k] + 0.4, dz) * (1 - smoothstep(24, 31, dz))
            : smoothstep(mL[k] - 0.4, mL[k] + 0.4, -dz) * (1 - smoothstep(19, 25, -dz));
          if (zg > 0 && gumK[k] > 0) {
            const qg = (bil(Ig, a00, a01, a10, a11, fu, fv) - fl) / (1 - fl);
            const hwg = Math.max(hw, qg > 0 ? hwMax * Math.pow(qg > 1 ? 1 : qg, gam) : 0, hwMax * 0.55 * Math.max(dentU[k], dentL[k]));
            if (hwg > 0.4) { let gw = (hwg + 1.7 - ad) / 1.25 + 0.5; gw = gw < 0 ? 0 : gw > 1 ? 1 : gw; B = gw * zg * gumK[k]; }
          }
          if (R + G + B <= 0) continue;
          const o = base + stride4 * kz;
          if (R > 0.25 || B > 0.5) { if (x < cmin[0]) cmin[0] = x; if (x > cmax[0]) cmax[0] = x; if (y < cmin[1]) cmin[1] = y; if (y > cmax[1]) cmax[1] = y; if (zc[kz] < cmin[2]) cmin[2] = zc[kz]; if (zc[kz] > cmax[2]) cmax[2] = zc[kz]; }
          rgba[o] = R >= 1 ? 255 : (R * 255 + 0.5) | 0;
          rgba[o + 1] = G >= 1 ? 255 : (G * 255 + 0.5) | 0;
          rgba[o + 2] = B >= 1 ? 255 : (B * 255 + 0.5) | 0;
          rgba[o + 3] = M >= 1 ? 255 : (M * 255 + 0.5) | 0;
        }
      }
    }
    blur3(rgba, nx, ny, nz);
    const slice = new Uint8Array(nvox);   // what the slice panes show: hard tissue plus a faint soft-tissue halo
    for (let i = 0, o = 0; i < nvox; i++, o += 4) { const r = rgba[o]; slice[i] = (r + (255 - r) * rgba[o + 2] * 0.17 / 255) | 0; }
    return { nx, ny, nz, vox, bmin, bmax, size: [bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]], Hmm, mmPerPx, rgba, slice, ms: performance.now() - t0, archW: P.archW, archD: P.archD, cmin: isFinite(cmin[0]) ? cmin : bmin.slice(), cmax: isFinite(cmax[0]) ? cmax : bmax.slice() };
  }
  /* ---- Teeth segmentation U-Net, plain JavaScript (mirrors teethseg.py, BatchNorm folded) ---- */
  function f16ToF32(u16) {
    const out = new Float32Array(u16.length);
    for (let i = 0; i < u16.length; i++) {
      const h = u16[i], s = h & 0x8000 ? -1 : 1, e = (h >> 10) & 31, f = h & 1023;
      out[i] = e === 0 ? s * f * 5.960464477539063e-8 : e === 31 ? (f ? NaN : s * Infinity) : s * Math.pow(2, e - 15) * (1 + f / 1024);
    }
    return out;
  }
  function makeNet(meta, u16) {
    const all = f16ToF32(u16);
    return { meta, layers: meta.layers.map(L => ({ cin: L.cin, cout: L.cout, k: L.k, w: all.subarray(L.w, L.w + L.cout * L.cin * L.k * L.k), b: all.subarray(L.b, L.b + L.cout) })) };
  }
  function conv3(x, h, w, L, relu) {
    // 4 output channels per pass so each input value is read once per tap
    const cin = L.cin, cout = L.cout, hw = h * w, out = new Float32Array(cout * hw), Wt = L.w, B = L.b;
    for (let o = 0; o < cout; o += 4) {
      const n = Math.min(4, cout - o), o0 = o * hw, o1 = o0 + hw, o2 = o1 + hw, o3 = o2 + hw;
      for (let q = 0; q < n; q++) out.fill(B[o + q], o0 + q * hw, o0 + (q + 1) * hw);
      for (let i = 0; i < cin; i++) {
        const xi = i * hw;
        for (let ky = 0; ky < 3; ky++) {
          const dy = ky - 1, y0 = dy < 0 ? 1 : 0, y1 = dy > 0 ? h - 1 : h;
          for (let kx = 0; kx < 3; kx++) {
            const t = ky * 3 + kx, dx = kx - 1, x0 = dx < 0 ? 1 : 0, x1 = dx > 0 ? w - 1 : w;
            const w0 = Wt[(o * cin + i) * 9 + t];
            const w1 = n > 1 ? Wt[((o + 1) * cin + i) * 9 + t] : 0;
            const w2 = n > 2 ? Wt[((o + 2) * cin + i) * 9 + t] : 0;
            const w3 = n > 3 ? Wt[((o + 3) * cin + i) * 9 + t] : 0;
            for (let y = y0; y < y1; y++) {
              const r = y * w, irow = xi + (y + dy) * w + dx;
              if (n === 4) {
                for (let xx = x0; xx < x1; xx++) { const v = x[irow + xx], k = r + xx; out[o0 + k] += w0 * v; out[o1 + k] += w1 * v; out[o2 + k] += w2 * v; out[o3 + k] += w3 * v; }
              } else {
                for (let xx = x0; xx < x1; xx++) { const v = x[irow + xx], k = r + xx; out[o0 + k] += w0 * v; if (n > 1) out[o1 + k] += w1 * v; if (n > 2) out[o2 + k] += w2 * v; }
              }
            }
          }
        }
      }
      if (relu) for (let k = o0; k < o0 + n * hw; k++) if (out[k] < 0) out[k] = 0;
    }
    return out;
  }
  function conv1(x, h, w, L) {
    const cin = L.cin, cout = L.cout, hw = h * w, out = new Float32Array(cout * hw);
    for (let o = 0; o < cout; o++) {
      const oo = o * hw; out.fill(L.b[o], oo, oo + hw);
      for (let i = 0; i < cin; i++) { const wv = L.w[o * cin + i], xi = i * hw; for (let k = 0; k < hw; k++) out[oo + k] += wv * x[xi + k]; }
    }
    return out;
  }
  function maxpool2(x, c, h, w) {
    const H = h >> 1, W = w >> 1, out = new Float32Array(c * H * W);
    for (let ch = 0; ch < c; ch++) for (let y = 0; y < H; y++) for (let xx = 0; xx < W; xx++) {
      const b = ch * h * w + 2 * y * w + 2 * xx;
      out[ch * H * W + y * W + xx] = Math.max(x[b], x[b + 1], x[b + w], x[b + w + 1]);
    }
    return out;
  }
  function up2cat(x, c, h, w, skip, cs) {
    const H = h * 2, W = w * 2, out = new Float32Array((c + cs) * H * W);
    for (let ch = 0; ch < c; ch++) for (let y = 0; y < H; y++) {
      const src = ch * h * w + (y >> 1) * w, dst = ch * H * W + y * W;
      for (let xx = 0; xx < W; xx++) out[dst + xx] = x[src + (xx >> 1)];
    }
    out.set(skip, c * H * W);
    return out;
  }
  function runUNet(net, inp, h, w) {
    const ch = net.meta.ch, nl = ch.length, off = net.meta.input_offset || 0;
    let li = 0, x = new Float32Array(inp.length), hh = h, ww = w, c = 1;
    for (let i = 0; i < inp.length; i++) x[i] = inp[i] + off;
    const skips = [];
    for (let i = 0; i < nl; i++) {
      x = conv3(x, hh, ww, net.layers[li++], true); x = conv3(x, hh, ww, net.layers[li++], true); c = ch[i];
      if (i < nl - 1) { skips.push({ x, c }); x = maxpool2(x, c, hh, ww); hh >>= 1; ww >>= 1; }
    }
    for (let i = nl - 2; i >= 0; i--) {
      const s = skips[i]; x = up2cat(x, c, hh, ww, s.x, s.c); hh <<= 1; ww <<= 1;
      x = conv3(x, hh, ww, net.layers[li++], true); x = conv3(x, hh, ww, net.layers[li++], true); c = ch[i];
    }
    const lg = conv1(x, hh, ww, net.layers[li++]), p = new Float32Array(hh * ww);
    for (let k = 0; k < p.length; k++) p[k] = 1 / (1 + Math.exp(-lg[k]));
    return p;
  }
  // Area-average resample (matches PIL's BOX filter closely enough for this model)
  function areaResample(src, W, H, w, h) {
    const tmp = new Float32Array(w * H), out = new Float32Array(w * h), sx = W / w, sy = H / h;
    for (let x = 0; x < w; x++) {
      const a = x * sx, b = a + sx, i0 = Math.floor(a), i1 = Math.min(W, Math.ceil(b));
      for (let y = 0; y < H; y++) {
        let s = 0; const r = y * W;
        for (let i = i0; i < i1; i++) s += src[r + i] * (Math.min(b, i + 1) - Math.max(a, i));
        tmp[y * w + x] = s / sx;
      }
    }
    for (let y = 0; y < h; y++) {
      const a = y * sy, b = a + sy, j0 = Math.floor(a), j1 = Math.min(H, Math.ceil(b));
      for (let x = 0; x < w; x++) {
        let s = 0;
        for (let j = j0; j < j1; j++) s += tmp[j * w + x] * (Math.min(b, j + 1) - Math.max(a, j));
        out[y * w + x] = s / sy;
      }
    }
    return out;
  }
  // norm: percentile-normalised OPG (W x H). Returns the teeth probability at network size.
  function segmentTeeth(net, norm, W, H) {
    const h = net.meta.in_h, w = net.meta.in_w;
    return { prob: runUNet(net, areaResample(norm, W, H, w, h), h, w), w, h };
  }

  /* ---- Surface extraction: naive surface nets on the RGBA volume, with gradient normals and baked AO ---- */
  const AO_D = new Float32Array([1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1].concat(
    [1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1].map(v => v / Math.sqrt(3))));
  const AO_R = [2, 4, 7], AO_W = [0.5, 0.3, 0.2];
  const CUBE_EDGES = [0, 1, 2, 3, 4, 5, 6, 7, 0, 2, 1, 3, 4, 6, 5, 7, 0, 4, 1, 5, 2, 6, 3, 7];

  // field: Int16Array (inside where > 0). Returns vertex positions in grid units and triangles.
  function surfaceNets(field, nx, ny, nz, lo, hi) {
    const sxy = nx * ny, layer = nx * ny;
    const buf = new Int32Array(2 * layer);
    let vcap = 1 << 16, fcap = 1 << 17, vc = 0, fc = 0;
    let V = new Float32Array(vcap * 3), F = new Uint32Array(fcap * 3);
    const g = new Float32Array(8);
    const x0 = Math.max(0, lo[0] - 1), y0 = Math.max(0, lo[1] - 1), z0 = Math.max(0, lo[2] - 1);
    const x1 = Math.min(nx - 1, hi[0] + 1), y1 = Math.min(ny - 1, hi[1] + 1), z1 = Math.min(nz - 1, hi[2] + 1);
    const addTri = (a, b, c) => {
      if (fc >= fcap) { fcap *= 2; const n = new Uint32Array(fcap * 3); n.set(F); F = n; }
      F[fc * 3] = a; F[fc * 3 + 1] = b; F[fc * 3 + 2] = c; fc++;
    };
    for (let z = z0; z < z1; z++) {
      const cur = (z & 1) * layer, prv = ((z + 1) & 1) * layer;
      for (let y = y0; y < y1; y++) {
        for (let x = x0; x < x1; x++) {
          const i = x + nx * y + sxy * z;
          const f0 = field[i], f1 = field[i + 1], f2 = field[i + nx], f3 = field[i + nx + 1];
          const f4 = field[i + sxy], f5 = field[i + sxy + 1], f6 = field[i + sxy + nx], f7 = field[i + sxy + nx + 1];
          const mask = (f0 > 0 ? 1 : 0) | (f1 > 0 ? 2 : 0) | (f2 > 0 ? 4 : 0) | (f3 > 0 ? 8 : 0) | (f4 > 0 ? 16 : 0) | (f5 > 0 ? 32 : 0) | (f6 > 0 ? 64 : 0) | (f7 > 0 ? 128 : 0);
          if (mask === 0 || mask === 255) continue;
          g[0] = f0; g[1] = f1; g[2] = f2; g[3] = f3; g[4] = f4; g[5] = f5; g[6] = f6; g[7] = f7;
          let px = 0, py = 0, pz = 0, n = 0;
          for (let e = 0; e < 24; e += 2) {
            const a = CUBE_EDGES[e], b = CUBE_EDGES[e + 1], ga = g[a], gb = g[b];
            if ((ga > 0) === (gb > 0)) continue;
            const t = ga / (ga - gb);
            px += (a & 1) + ((b & 1) - (a & 1)) * t;
            py += ((a >> 1) & 1) + (((b >> 1) & 1) - ((a >> 1) & 1)) * t;
            pz += ((a >> 2) & 1) + (((b >> 2) & 1) - ((a >> 2) & 1)) * t;
            n++;
          }
          if (vc >= vcap) { vcap *= 2; const nv = new Float32Array(vcap * 3); nv.set(V); V = nv; }
          V[vc * 3] = x + px / n; V[vc * 3 + 1] = y + py / n; V[vc * 3 + 2] = z + pz / n;
          const vi = vc++;
          buf[cur + x + nx * y] = vi;
          // one quad per crossing edge leaving corner 0, built from the four cells that share it
          if ((mask & 1) !== ((mask >> 1) & 1) && y > y0 && z > z0) {          // x-edge
            const a = vi, b = buf[cur + x + nx * (y - 1)], c = buf[prv + x + nx * (y - 1)], d = buf[prv + x + nx * y];
            if (mask & 1) { addTri(a, b, c); addTri(a, c, d); } else { addTri(a, d, c); addTri(a, c, b); }
          }
          if ((mask & 1) !== ((mask >> 2) & 1) && x > x0 && z > z0) {          // y-edge
            const a = vi, b = buf[prv + x + nx * y], c = buf[prv + (x - 1) + nx * y], d = buf[cur + (x - 1) + nx * y];
            if (mask & 1) { addTri(a, b, c); addTri(a, c, d); } else { addTri(a, d, c); addTri(a, c, b); }
          }
          if ((mask & 1) !== ((mask >> 4) & 1) && x > x0 && y > y0) {          // z-edge
            const a = vi, b = buf[cur + (x - 1) + nx * y], c = buf[cur + (x - 1) + nx * (y - 1)], d = buf[cur + x + nx * (y - 1)];
            if (mask & 1) { addTri(a, b, c); addTri(a, c, d); } else { addTri(a, d, c); addTri(a, c, b); }
          }
        }
      }
    }
    return { V: V.subarray(0, vc * 3), F: F.subarray(0, fc * 3), vc, fc };
  }
  function sampleField(f, nx, ny, nz, x, y, z) {
    // trilinear on the grid (grid coordinates), clamped
    x = x < 0 ? 0 : x > nx - 1.001 ? nx - 1.001 : x; y = y < 0 ? 0 : y > ny - 1.001 ? ny - 1.001 : y; z = z < 0 ? 0 : z > nz - 1.001 ? nz - 1.001 : z;
    const x0 = x | 0, y0 = y | 0, z0 = z | 0, fx = x - x0, fy = y - y0, fz = z - z0, sxy = nx * ny, i = x0 + nx * y0 + sxy * z0;
    const c00 = f[i] + (f[i + 1] - f[i]) * fx, c10 = f[i + nx] + (f[i + nx + 1] - f[i + nx]) * fx;
    const c01 = f[i + sxy] + (f[i + sxy + 1] - f[i + sxy]) * fx, c11 = f[i + sxy + nx] + (f[i + sxy + nx + 1] - f[i + sxy + nx]) * fx;
    const c0 = c00 + (c10 - c00) * fy, c1 = c01 + (c11 - c01) * fy;
    return c0 + (c1 - c0) * fz;
  }
  // cheap 3D value noise for a little surface grain (baked per vertex)
  function hash3(x, y, z) { let h = Math.imul(x, 374761393) ^ Math.imul(y, 668265263) ^ Math.imul(z, 2147483647); h = Math.imul(h ^ (h >>> 13), 1274126177); return ((h ^ (h >>> 16)) >>> 0) / 4294967296; }
  function vnoise(x, y, z) {
    const xi = Math.floor(x), yi = Math.floor(y), zi = Math.floor(z);
    let fx = x - xi, fy = y - yi, fz = z - zi; fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy); fz = fz * fz * (3 - 2 * fz);
    const L = (a, b, t) => a + (b - a) * t;
    return L(L(L(hash3(xi, yi, zi), hash3(xi + 1, yi, zi), fx), L(hash3(xi, yi + 1, zi), hash3(xi + 1, yi + 1, zi), fx), fy),
      L(L(hash3(xi, yi, zi + 1), hash3(xi + 1, yi, zi + 1), fx), L(hash3(xi, yi + 1, zi + 1), hash3(xi + 1, yi + 1, zi + 1), fx), fy), fz);
  }
  // kind: 'hard' (attr.x = toothness), 'teeth' (attr.x = crown, i.e. not buried in bone), 'gum'
  function finishMesh(net, field, rgba, occ, occBits, kind, T, nx, ny, nz, bmin, vox) {
    const { V, F, vc } = net, sxy = nx * ny;
    const AO_O = new Int32Array(14 * 3 * 3), AO_L = new Int32Array(14 * 3);
    for (let k = 0; k < 14; k++) for (let r = 0; r < 3; r++) {
      const q = (k * 3 + r) * 3;
      for (let a = 0; a < 3; a++) AO_O[q + a] = Math.round(AO_D[k * 3 + a] * AO_R[r]);
      AO_L[k * 3 + r] = AO_O[q] + nx * AO_O[q + 1] + sxy * AO_O[q + 2];
    }
    const pos = new Float32Array(vc * 3), nrm = new Int8Array(vc * 4), attr = new Uint8Array(vc * 4);
    const at = (x, y, z) => {
      const ix = x < 0 ? 0 : x > nx - 1 ? nx - 1 : Math.round(x), iy = y < 0 ? 0 : y > ny - 1 ? ny - 1 : Math.round(y), iz = z < 0 ? 0 : z > nz - 1 ? nz - 1 : Math.round(z);
      return (ix + nx * iy + sxy * iz) * 4;
    };
    for (let v = 0; v < vc; v++) {
      const gx = V[v * 3], gy = V[v * 3 + 1], gz = V[v * 3 + 2];
      const wx = bmin[0] + (gx + 0.5) * vox, wy = bmin[1] + (gy + 0.5) * vox, wz = bmin[2] + (gz + 0.5) * vox;
      pos[v * 3] = wx; pos[v * 3 + 1] = wy; pos[v * 3 + 2] = wz;
      let dx = sampleField(field, nx, ny, nz, gx - 1.2, gy, gz) - sampleField(field, nx, ny, nz, gx + 1.2, gy, gz);
      let dy = sampleField(field, nx, ny, nz, gx, gy - 1.2, gz) - sampleField(field, nx, ny, nz, gx, gy + 1.2, gz);
      let dz = sampleField(field, nx, ny, nz, gx, gy, gz - 1.2) - sampleField(field, nx, ny, nz, gx, gy, gz + 1.2);
      const l = Math.hypot(dx, dy, dz) || 1; dx /= l; dy /= l; dz /= l;          // outward normal (field decreases outward)
      nrm[v * 4] = Math.round(dx * 127); nrm[v * 4 + 1] = Math.round(dy * 127); nrm[v * 4 + 2] = Math.round(dz * 127);
      const o = at(gx - dx * 0.8, gy - dy * 0.8, gz - dz * 0.8);            // material just inside the surface
      if (kind === 'hard') attr[v * 4] = Math.max(rgba[o + 1], rgba[o + 3]);
      else if (kind === 'teeth') { const q = at(gx + dx * 1.6, gy + dy * 1.6, gz + dz * 1.6); attr[v * 4] = rgba[q] >= T && rgba[q + 1] < 128 ? 40 : 255; }
      else attr[v * 4] = 0;
      attr[v * 4 + 1] = kind === 'gum' ? 0 : rgba[o + 3];
      // ambient occlusion: how much solid tissue surrounds the point on its outer side (cosine-weighted)
      const bx = (gx + 0.5) | 0, by = (gy + 0.5) | 0, bz = (gz + 0.5) | 0, bi = bx + nx * by + sxy * bz;
      let occl = 0, tot = 0;
      for (let k = 0; k < 14; k++) {
        const c = AO_D[k * 3] * dx + AO_D[k * 3 + 1] * dy + AO_D[k * 3 + 2] * dz; if (c < 0.2) continue;
        tot += c * 1.4;
        for (let r = 0; r < 3; r++) {
          const q = (k * 3 + r) * 3, sx = bx + AO_O[q], sy = by + AO_O[q + 1], sz = bz + AO_O[q + 2];
          if (sx < 0 || sy < 0 || sz < 0 || sx >= nx || sy >= ny || sz >= nz) continue;
          if (occ[bi + AO_L[k * 3 + r]] & occBits) occl += c * AO_W[r];
        }
      }
      attr[v * 4 + 2] = Math.round(255 * Math.max(0.08, 1 - 1.6 * occl / (tot || 1)));
      const nz1 = vnoise(wx * 0.55, wy * 0.55, wz * 0.55), nz2 = vnoise(wx * 1.7 + 17, wy * 1.7, wz * 1.7);
      attr[v * 4 + 3] = Math.round(255 * (0.65 * nz1 + 0.35 * nz2));
    }
    return { pos, nrm, attr, idx: F.slice(), vc, fc: F.length / 3 };
  }
  function buildMeshes(rgba, nx, ny, nz, bmin, vox, thr, cmin, cmax) {
    const t0 = performance.now(), n = nx * ny * nz, T = Math.round(thr * 255), sxy = nx * ny;
    const lo = [0, 1, 2].map(a => Math.max(0, Math.floor((cmin[a] - bmin[a]) / vox) - 3));
    const hi = [nx, ny, nz].map((d, a) => Math.min(d - 1, Math.ceil((cmax[a] - bmin[a]) / vox) + 3));
    const field = new Int16Array(n), occ = new Uint8Array(n);
    const X0 = Math.max(0, lo[0] - 2), Y0 = Math.max(0, lo[1] - 2), Z0 = Math.max(0, lo[2] - 2);
    const X1 = Math.min(nx - 1, hi[0] + 2), Y1 = Math.min(ny - 1, hi[1] + 2), Z1 = Math.min(nz - 1, hi[2] + 2);
    for (let z = Z0; z <= Z1; z++) for (let y = Y0; y <= Y1; y++) {
      for (let i = X0 + nx * y + sxy * z, o = i * 4, e = X1 + nx * y + sxy * z; i <= e; i++, o += 4) {
        const hard = rgba[o] >= T, tooth = hard && (rgba[o + 1] >= 128 || rgba[o + 3] >= 128);
        occ[i] = (hard ? 1 : 0) | (rgba[o + 2] >= 128 ? 2 : 0) | (tooth ? 4 : 0);
      }
    }
    const out = { thr };
    const bits = { hard: 1, gum: 3, teeth: 4 };
    for (const kind of ['hard', 'teeth', 'gum']) {
      const K = kind === 'hard' ? 0 : kind === 'gum' ? 2 : 1;
      field.fill(-1);
      for (let z = Z0; z <= Z1; z++) for (let y = Y0; y <= Y1; y++) {
        const edgeYZ = z === 0 || y === 0 || z === nz - 1 || y === ny - 1;
        for (let x = X0, i = X0 + nx * y + sxy * z, o = i * 4; x <= X1; x++, i++, o += 4) {
          if (edgeYZ || x === 0 || x === nx - 1) continue;           // closes the surface at the grid edge
          if (K === 0) field[i] = rgba[o] - T;
          else if (K === 2) field[i] = rgba[o + 2] - 128;
          else { const r = rgba[o] - T, t = (rgba[o + 1] > rgba[o + 3] ? rgba[o + 1] : rgba[o + 3]) - 128; field[i] = r < t ? r : t; }
        }
      }
      out[kind] = finishMesh(surfaceNets(field, nx, ny, nz, lo, hi), field, rgba, occ, bits[kind], kind, T, nx, ny, nz, bmin, vox);
    }
    out.ms = performance.now() - t0;
    return out;
  }

  function teethToMap(t, W, H) {
    const out = new Float32Array(W * H), sx = t.w / W, sy = t.h / H, P = t.prob;
    for (let y = 0; y < H; y++) {
      const fy = (y + 0.5) * sy - 0.5, y0 = clamp(Math.floor(fy), 0, t.h - 1), y1 = Math.min(y0 + 1, t.h - 1), wy = clamp(fy - y0, 0, 1);
      for (let x = 0; x < W; x++) {
        const fx = (x + 0.5) * sx - 0.5, x0 = clamp(Math.floor(fx), 0, t.w - 1), x1 = Math.min(x0 + 1, t.w - 1), wx = clamp(fx - x0, 0, 1);
        out[y * W + x] = (P[y0 * t.w + x0] * (1 - wx) + P[y0 * t.w + x1] * wx) * (1 - wy) + (P[y1 * t.w + x0] * (1 - wx) + P[y1 * t.w + x1] * wx) * wy;
      }
    }
    return out;
  }
  return { clamp, lerp, smoothstep, UPPER, LOWER, REF_L, thickAt, buildArch, nearestIdx, idxAtS, tableAt, occFracAt, gaussBlur, blurXY, percentileNormalize, buildVolume, makeNet: typeof makeNet === 'function' ? makeNet : null, segmentTeeth: typeof segmentTeeth === 'function' ? segmentTeeth : null, buildMeshes: typeof buildMeshes === 'function' ? buildMeshes : null };
}

function workerMain() {
  const SH = shared(); let src = null, net = null, last = null; const cache = {};
  const bufs = m => { const out = []; for (const k of ['hard', 'teeth', 'gum']) { const g = m[k]; out.push(g.pos.buffer, g.nrm.buffer, g.attr.buffer, g.idx.buffer); } return out; };
  const fail = (id, err) => self.postMessage({ type: 'error', id, msg: String((err && err.message) || err) });
  self.onmessage = e => {
    const m = e.data;
    if (m.type === 'net') { net = SH.makeNet(m.meta, m.u16); return; }
    if (m.type === 'segment') {                    // runs in its own worker so the 3D stays responsive
      if (!net) { self.postMessage({ type: 'teeth', sid: m.sid, error: 'no model' }); return; }
      try {
        const t0 = performance.now(), t = SH.segmentTeeth(net, m.norm, m.w, m.h);
        self.postMessage({ type: 'teeth', sid: m.sid, prob: t.prob, w: t.w, h: t.h, ms: performance.now() - t0 }, [t.prob.buffer]);
      } catch (err) { self.postMessage({ type: 'teeth', sid: m.sid, error: String((err && err.message) || err) }); }
      return;
    }
    if (m.type === 'source') { src = m; src.teeth = null; cache.key = null; last = null; return; }
    if (m.type === 'teeth') { if (src && src.sid === m.sid) src.teeth = m.prob ? { prob: m.prob, w: m.w, h: m.h } : null; return; }
    if (m.type === 'build') {
      try {
        const r = SH.buildVolume(src.norm, src.raw, src.w, src.h, m.params, cache, m.params.teethSrc === 'model' ? src.teeth : null);
        r.mesh = SH.buildMeshes(r.rgba, r.nx, r.ny, r.nz, r.bmin, r.vox, m.thr, r.cmin, r.cmax);
        last = { rgba: r.rgba.slice(), nx: r.nx, ny: r.ny, nz: r.nz, bmin: r.bmin, vox: r.vox, cmin: r.cmin, cmax: r.cmax, ms: r.ms };
        self.postMessage({ type: 'built', id: m.id, r }, [r.rgba.buffer, r.slice.buffer].concat(bufs(r.mesh)));
      } catch (err) { fail(m.id, err); }
      return;
    }
    if (m.type === 'mesh') {
      if (!last || last.ms !== m.vms) return;
      try { const mesh = SH.buildMeshes(last.rgba, last.nx, last.ny, last.nz, last.bmin, last.vox, m.thr, last.cmin, last.cmax); self.postMessage({ type: 'meshed', id: m.id, vms: m.vms, mesh }, bufs(mesh)); }
      catch (err) { self.postMessage({ type: 'meshed', id: m.id, error: String(err) }); }
    }
  };
}

const SH = shared();
const { clamp, lerp, smoothstep, UPPER, LOWER, REF_L, thickAt, nearestIdx, idxAtS, occFracAt } = SH;
const $ = (s, r = document) => r.querySelector(s);
const V3 = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
  sub: (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]],
  mul: (a, s) => [a[0] * s, a[1] * s, a[2] * s],
  dot: (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2],
  cross: (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]],
  norm: a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; },
  len: a => Math.hypot(a[0], a[1], a[2]),
};
const mm = v => (v < 0 ? '−' : '') + Math.abs(v).toFixed(1);
const TOOTH_NAMES = ['central incisor', 'lateral incisor', 'canine', 'first premolar', 'second premolar', 'first molar', 'second molar', 'third molar'];

/* =====================================================================
   Synthetic demo OPG (drawn procedurally, no patient data)
   ===================================================================== */
function mulberry32(a) { return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
const DEMO = { W: 1400, H: 700 };
DEMO.ppm = DEMO.W / (2 * REF_L);
DEMO.occY = x => { const n = (x - DEMO.W / 2) / (DEMO.W / 2); return DEMO.H * 0.515 - DEMO.H * 0.075 * n * n; };
DEMO.occFracAtS = s => DEMO.occY(DEMO.W / 2 + s * DEMO.ppm) / DEMO.H;

function makeDemoOPG() {
  const { W, H, ppm } = DEMO, cx = W / 2, hw = W / 2, gaussBlur = SH.gaussBlur;
  const occY = DEMO.occY;
  const occSlope = x => -2 * H * 0.075 * ((x - cx) / hw) / hw;
  const cv = document.createElement('canvas'); cv.width = W; cv.height = H;
  const g = cv.getContext('2d', { willReadFrequently: true });
  const mask = fn => {
    g.setTransform(1, 0, 0, 1, 0, 0); g.fillStyle = '#000'; g.fillRect(0, 0, W, H);
    g.fillStyle = '#fff'; g.strokeStyle = '#fff'; g.lineCap = 'round'; g.lineJoin = 'round';
    fn(g);
    const d = g.getImageData(0, 0, W, H).data, m = new Float32Array(W * H);
    for (let i = 0, j = 0; i < m.length; i++, j += 4) m[i] = d[j] / 255;
    return m;
  };
  const P = (n, y) => [cx + n * hw, y];
  const mirrorLoop = right => right.map(([n, y]) => P(n, y)).concat(right.slice().reverse().map(([n, y]) => P(-n, y)));
  const poly = pts => { g.beginPath(); pts.forEach((p, i) => i ? g.lineTo(p[0], p[1]) : g.moveTo(p[0], p[1])); g.closePath(); };
  const oy = n => occY(cx + n * hw);
  const mand = [];
  for (let n = 0; n <= 0.741; n += 0.037) mand.push([n, oy(n) + lerp(88, 74, smoothstep(0.12, 0.45, n))]);
  mand.push([0.775, oy(0.775) + 48], [0.795, oy(0.795) - 5], [0.803, 235], [0.8, 178], [0.806, 160], [0.822, 168], [0.845, 200], [0.872, 214], [0.895, 196], [0.905, 160], [0.912, 126], [0.93, 106], [0.952, 110], [0.962, 132], [0.957, 168], [0.964, 250], [0.968, 370], [0.96, 470], [0.94, 552], [0.905, 596], [0.83, 620], [0.66, 640], [0.45, 655], [0.22, 665], [0, 669]);
  const maxl = [];
  for (let n = 0; n <= 0.741; n += 0.037) maxl.push([n, oy(n) - lerp(98, 76, smoothstep(0.12, 0.45, n))]);
  maxl.push([0.772, oy(0.772) - 50], [0.79, 255], [0.785, 190], [0.75, 120], [0.62, 75], [0.42, 58], [0.2, 62], [0, 66]);
  const mMand = mask(() => { poly(mirrorLoop(mand)); g.fill(); });
  const mMax = mask(() => { poly(mirrorLoop(maxl)); g.fill(); });
  const mSoft = gaussBlur(mask(() => { g.beginPath(); g.ellipse(cx, H * 0.56, hw * 1.02, H * 0.62, 0, 0, 7); g.fill(); }), W, H, 40);
  const mSpine = gaussBlur(mask(() => g.fillRect(cx - 60, H * 0.25, 120, H)), W, H, 30);
  const sinusPath = s => { g.beginPath(); g.ellipse(cx + s * 0.44 * hw, 168, 0.18 * hw, 70, s * 0.08, 0, 7); };
  const mSinus = gaussBlur(mask(() => { for (const s of [-1, 1]) { sinusPath(s); g.fill(); } }), W, H, 3);
  const mSinRim = gaussBlur(mask(() => { g.lineWidth = 4; for (const s of [-1, 1]) { sinusPath(s); g.stroke(); } }), W, H, 1.2);
  const mNasal = gaussBlur(mask(() => { g.beginPath(); g.ellipse(cx, 120, 0.15 * hw, 88, 0, 0, 7); g.fill(); }), W, H, 6);
  const mSeptum = gaussBlur(mask(() => { g.lineWidth = 5; g.beginPath(); g.moveTo(cx, 40); g.lineTo(cx, 200); g.stroke(); }), W, H, 1.5);
  const mPalate = gaussBlur(mask(() => { g.lineWidth = 13; g.beginPath(); g.moveTo(cx - 0.5 * hw, 218); g.quadraticCurveTo(cx, 196, cx + 0.5 * hw, 218); g.stroke(); }), W, H, 3);
  const mZyg = gaussBlur(mask(() => { g.lineWidth = 16; for (const s of [-1, 1]) { g.beginPath(); g.moveTo(cx + s * 0.6 * hw, 262); g.quadraticCurveTo(cx + s * 0.66 * hw, 190, cx + s * 0.78 * hw, 128); g.stroke(); } }), W, H, 5);
  const canalPath = s => { g.beginPath(); g.moveTo(cx + s * 0.885 * hw, 330); g.quadraticCurveTo(cx + s * 0.84 * hw, 540, cx + s * 0.62 * hw, 566); g.quadraticCurveTo(cx + s * 0.45 * hw, 578, cx + s * 0.33 * hw, 572); };
  const mCanal = gaussBlur(mask(() => { g.lineWidth = 9; for (const s of [-1, 1]) { canalPath(s); g.stroke(); } }), W, H, 1.5);
  const mCanalW = gaussBlur(mask(() => { g.lineWidth = 15; for (const s of [-1, 1]) { canalPath(s); g.stroke(); } }), W, H, 1.2);
  const mandBlur = gaussBlur(mMand, W, H, 5);
  const rnd = mulberry32(7);
  const tr = new Float32Array(W * H); for (let i = 0; i < tr.length; i++) tr[i] = rnd() * 2 - 1;
  const trab = gaussBlur(tr, W, H, 1.6);
  const img = new Float32Array(W * H);
  for (let i = 0; i < img.length; i++) {
    const y = (i / W) | 0;
    const cortex = clamp((mMand[i] - mandBlur[i]) * 2.2, 0, 1);
    const mx = mMax[i] * smoothstep(70, 175, y);
    let v = 0.07 + 0.12 * mSoft[i] + 0.06 * mSpine[i];
    v += mx * 0.26 * (1 - 0.62 * mSinus[i]) * (1 - 0.5 * mNasal[i]) + 0.1 * mSinRim[i] * mMax[i] + 0.16 * mPalate[i] + 0.06 * mSeptum[i] + 0.1 * mZyg[i] * mMax[i];
    v += mMand[i] * 0.27 + 0.26 * cortex - 0.1 * mCanal[i] * mMand[i] + 0.06 * Math.max(0, mCanalW[i] - mCanal[i]) * mMand[i];
    v += (mMand[i] + mx) * trab[i] * 0.22;
    img[i] = v;
  }
  const special = {
    16: { amalgam: 'mod' }, 36: { amalgam: 'occ' }, 46: { rct: true, metalCrown: true },
    25: { implant: true, metalCrown: true }, 38: { tilt: -0.75, drop: 5, shift: 1.5 }, 28: { missing: true },
  };
  for (const up of [true, false]) for (const sg of [-1, 1]) {
    const T = up ? UPPER : LOWER;
    for (let n = 0; n < 8; n++) {
      const q = up ? (sg < 0 ? 1 : 2) : (sg < 0 ? 4 : 3);
      const sp = special[q * 10 + n + 1] || {};
      if (sp.missing) continue;
      drawTooth(Object.assign({ up, s: sg * (T[n].c + (sp.shift || 0)) }, T[n], sp));
    }
  }
  function drawTooth(t) {
    const x0 = cx + t.s * ppm;
    const ang = Math.atan(occSlope(x0)) + (t.tilt || 0);
    const ca = Math.cos(ang), sa = Math.sin(ang);
    const ap = t.up ? [sa, -ca] : [-sa, ca], md = [ca, sa];
    const chp = t.ch * ppm, cwp = (t.w - 0.5) * ppm, rlp = t.rl * ppm;
    const off = 3 + (t.drop || 0) * ppm + chp;
    const c = [x0 + ap[0] * off, occY(x0) + ap[1] * off];
    const roots = [];
    if (t.implant) { /* body drawn below */ }
    else if (t.r === 1) roots.push({ o0: 0, o1: 0, hb: cwp * 0.33, len: rlp });
    else if (t.r === 2 && t.w < 8) roots.push({ o0: -cwp * 0.12, o1: -cwp * 0.15, hb: cwp * 0.19, len: rlp }, { o0: cwp * 0.12, o1: cwp * 0.15, hb: cwp * 0.19, len: rlp * 0.96 });
    else {
      roots.push({ o0: -cwp * 0.24, o1: -cwp * 0.3, hb: cwp * 0.2, len: rlp }, { o0: cwp * 0.24, o1: cwp * 0.27, hb: cwp * 0.2, len: rlp * 0.95 });
      if (t.r === 3) roots.push({ o0: 0, o1: cwp * 0.02, hb: cwp * 0.17, len: rlp * 1.05, faint: 0.55 });
    }
    const reach = chp + rlp * 1.15 + 14;
    const xa = Math.max(0, Math.floor(c[0] - reach)), xb = Math.min(W - 1, Math.ceil(c[0] + reach));
    const ya = Math.max(0, Math.floor(c[1] - reach)), yb = Math.min(H - 1, Math.ceil(c[1] + reach));
    const cr = Math.min(cwp, chp) * 0.5;
    for (let y = ya; y <= yb; y++) for (let x = xa; x <= xb; x++) {
      const dx = x - c[0], dy = y - c[1];
      const lx = dx * md[0] + dy * md[1], ly = dx * ap[0] + dy * ap[1];
      let crown = 0, inner = 0, pulpC = 0;
      if (ly > -chp - 2 && ly < chp * 0.3) {
        const qx = Math.abs(lx) / (cwp / 2), qy = Math.abs(ly + chp * 0.45) / (chp * 0.55);
        crown = clamp((1 - Math.cbrt(qx * qx * qx + qy * qy * qy)) * cr + 0.5, 0, 1);
        const qx2 = Math.abs(lx) / (cwp / 2 - 0.9 * ppm), qy2 = Math.abs(ly + chp * 0.45) / (chp * 0.55 - 0.9 * ppm);
        inner = clamp((1 - Math.cbrt(qx2 * qx2 * qx2 + qy2 * qy2 * qy2)) * cr + 0.5, 0, 1);
        pulpC = clamp((1 - Math.hypot(lx / (cwp * 0.2), (ly + chp * 0.18) / (chp * 0.26))) * 6, 0, 1);
      }
      let root = 0, canal = 0, pdl = 0, lam = 0, implant = 0;
      for (const R of roots) {
        const tt = ly / R.len; if (tt < -0.12 || tt > 1.14) continue;
        const cxr = lerp(R.o0, R.o1, clamp(tt, 0, 1));
        const hwr = R.hb * (1 - 0.72 * Math.pow(clamp(tt, 0, 1), 1.2));
        const lat = Math.abs(lx - cxr);
        const dist = tt <= 1 ? lat - hwr : Math.hypot(lat, (tt - 1) * R.len) - R.hb * 0.28;
        const f = R.faint || 1;
        root = Math.max(root, clamp(0.5 - dist, 0, 1) * f);
        pdl = Math.max(pdl, clamp(0.5 - (dist - 2.2), 0, 1));
        lam = Math.max(lam, clamp(0.5 - (dist - 4.5), 0, 1));
        if (tt > 0 && tt < 0.96) { const cw = Math.max(0.9, (0.45 * (1 - tt) + 0.12) * ppm); canal = Math.max(canal, clamp(0.5 - (lat - cw), 0, 1) * f); }
      }
      if (t.implant) {
        const lm = ly / ppm;
        if (lm > -1.6 && lm < 11.9) {
          let r;
          if (lm < 0.4) r = 1.5;
          else {
            const tt = (lm - 0.4) / 11;
            r = 2.05 - 0.45 * tt * tt + 0.28 * (0.5 + 0.5 * Math.sin(lm * 2 * Math.PI / 0.85));
            if (lm > 11) r *= Math.sqrt(Math.max(0, 1 - ((lm - 11) / 0.9) ** 2));
          }
          implant = clamp(0.5 - (Math.abs(lx) / ppm - r) * ppm, 0, 1);
        }
      }
      if (crown <= 0 && root <= 0 && lam <= 0 && implant <= 0) continue;
      const tooth = Math.max(crown, root);
      const occl = clamp((-ly - chp * 0.35) / (chp * 0.65), 0, 1);
      let v = Math.max(crown * 0.36, root * 0.28) + crown * (0.08 + 0.12 * occl) + (crown - inner) * 0.08;
      const pulp = Math.max(pulpC * crown, canal) * tooth;
      if (t.implant) v += implant * 0.72;
      else if (t.rct) v += Math.max(canal * 0.7, pulpC * crown * 0.6);
      else v -= pulp * 0.2;
      if (t.metalCrown) v += crown * 0.5;
      if (t.amalgam) {
        const ax = lx / (cwp * (t.amalgam === 'mod' ? 0.47 : 0.24)), ay = (ly + chp * 0.8) / (chp * 0.2);
        v += clamp((1 - Math.hypot(ax, ay)) * 8, 0, 1) * crown * 0.6;
      }
      if (!t.implant) v += (-0.06 * Math.max(0, pdl - root) + 0.05 * Math.max(0, lam - pdl)) * (1 - crown);
      img[y * W + x] += v;
    }
  }
  for (let i = 0; i < img.length; i++) img[i] += (rnd() + rnd() + rnd() - 1.5) * 0.025;
  const out = gaussBlur(img, W, H, 0.9);
  for (let i = 0; i < out.length; i++) {
    const n = ((i % W) - cx) / hw;
    out[i] = Math.pow(clamp(out[i] * (0.93 + 0.07 * (1 - n ** 8)), 0, 1), 0.95);
  }
  return out;
}

/* =====================================================================
   State
   ===================================================================== */
/* Trained teeth-segmentation weights (fp16): either embedded as base64 or fetched from meta.src. */
const NETMETA = (() => { try { const m = document.getElementById('teethnet-meta'); return m && SH.makeNet ? JSON.parse(m.textContent) : null; } catch (e) { return null; } })();
let NETU16 = null;
const netReady = (async () => {
  if (!NETMETA) return false;
  try {
    const w = document.getElementById('teethnet-w');
    if (w && w.textContent.trim()) {
      const bin = atob(w.textContent.trim()), u8 = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
      NETU16 = new Uint16Array(u8.buffer);
    } else if (NETMETA.src) {
      const r = await fetch(NETMETA.src);
      if (!r.ok) throw new Error('HTTP ' + r.status);
      NETU16 = new Uint16Array(await r.arrayBuffer());
    } else return false;
    return true;
  } catch (e) { console.warn('Teeth model unavailable', e); return false; }
})();
let mainNet = null;
const getMainNet = () => mainNet || (NETU16 ? (mainNet = SH.makeNet(NETMETA, NETU16)) : null);
/* Phones and small or low-memory devices get a lighter model and fewer pixels. */
const LITE = /[?&]lite\b/.test(location.search) || (!/[?&]full\b/.test(location.search) && (
  (window.matchMedia && matchMedia('(pointer: coarse)').matches && Math.min(screen.width, screen.height) < 700) ||
  (navigator.deviceMemory && navigator.deviceMemory <= 2) || (navigator.hardwareConcurrency && navigator.hardwareConcurrency <= 2)));
const DEFAULT_PARAMS = { mid: 0.5, span: 1, archW: 1, archD: 1, vscale: 1, occl: 0.515, curve: 0.075, thick: 1, floor: 0.28, shell: 0.5, smooth: 1, round: 0.85, res: LITE ? 'lite' : 'std', invert: false, teethSrc: NETMETA ? 'model' : 'contrast' };
const PRESETS = {
  anatomy: { mode: 'surface', gums: true, boneA: 1 },
  bone: { mode: 'surface', gums: false, boneA: 1 },
  roots: { mode: 'surface', gums: false, boneA: 0.16 },
  xray: { mode: 'xray' },
};
const state = {
  src: null,
  params: { ...DEFAULT_PARAMS },
  render: { preset: 'anatomy', mode: 'surface', gums: true, gumA: 1, boneA: 1, thr: 0.33, clip: 'none', planes: false, spin: false },
  mpr: { win: 255, lev: 128, wpreset: 'full', obl: 'arch', yaw: 35, pitch: 20 },
  arch: null, vol: null, mesh: null,
  P: [0, 0, 0],
  cam: { yaw: -0.5, pitch: 0.22, zoom: 1, pan: [0, 0, 0] },
  moving: false,
  measure: false, measures: [],
  ui: { showTeeth: true },
};
const PANE_KEYS = ['axial', 'coronal', 'sagittal', 'oblique'];

/* =====================================================================
   Build pipeline (Web Worker with a same-thread fallback)
   ===================================================================== */
let worker = null, segWorker = null, buildSeq = 0, pending = null, current = null;
const mainCache = {};
function initWorker() {
  let url = null;
  try {
    url = URL.createObjectURL(new Blob([`const shared = ${shared.toString()};\n(${workerMain.toString()})();`], { type: 'text/javascript' }));
    worker = new Worker(url);
    worker.onmessage = e => {
      const m = e.data;
      if (m.type === 'meshed') { onMeshed(m); return; }
      if (!current || m.id !== current.id) return;
      if (m.type === 'built') finishBuild(m.r, current);
      else { console.error(m.msg); worker = null; runJob(current); }
    };
    worker.onerror = e => {
      console.warn('Worker unavailable, working on the main thread.', e.message || ''); e.preventDefault && e.preventDefault(); worker = null;
      if (current) runJob(current);
    };
  } catch (e) { worker = null; }
  if (!NETMETA || !url) return;
  try {
    segWorker = new Worker(url);
    netReady.then(ok => { if (ok && segWorker) segWorker.postMessage({ type: 'net', meta: NETMETA, u16: NETU16.slice() }); });
    segWorker.onmessage = e => { if (e.data.type === 'teeth') onTeeth(e.data); };
    segWorker.onerror = e => {
      e.preventDefault && e.preventDefault(); segWorker = null;
      const S = state.src; if (S && S.teethPending) netReady.then(() => setTimeout(segmentOnMain, 30));
    };
  } catch (e) { segWorker = null; }
}
function wantModel() { return !!NETMETA && state.params.teethSrc === 'model'; }
/* Hand the (possibly inverted) image to the build worker. Teeth from the model follow separately. */
function sendSource() {
  const S = state.src; if (!S || !worker) return;
  worker.postMessage({ type: 'source', sid: S.sid, norm: S.norm, raw: S.prep, w: S.w, h: S.h });
  if (S.teeth) worker.postMessage({ type: 'teeth', sid: S.sid, prob: S.teeth.prob, w: S.teeth.w, h: S.teeth.h });
}
function requestSegment() {
  const S = state.src; if (!S) return;
  S.teethPending = true; S.teeth = null; setStage();
  netReady.then(ok => {
    if (state.src !== S) return;
    if (!ok) { onTeeth({ sid: S.sid, error: 'no model' }); return; }
    if (segWorker) segWorker.postMessage({ type: 'segment', sid: S.sid, norm: S.norm, w: S.w, h: S.h });
    else setTimeout(segmentOnMain, 30);
  });
}
function segmentOnMain() {
  const S = state.src, net = getMainNet(); if (!S) return;
  if (!net) { onTeeth({ sid: S.sid, error: 'no model' }); return; }
  const t0 = performance.now();
  try { const t = SH.segmentTeeth(net, S.norm, S.w, S.h); onTeeth({ sid: S.sid, prob: t.prob, w: t.w, h: t.h, ms: performance.now() - t0 }); }
  catch (e) { console.error(e); onTeeth({ sid: S.sid, error: String(e) }); }
}
/* Teeth probability arrived: tint the thumbnail, pass it to the builder, then run whatever was waiting on it. */
function onTeeth(m) {
  const S = state.src; if (!S || m.sid !== S.sid) return;
  S.teethPending = false;
  if (m.error) { S.teeth = null; state.params.teethSrc = 'contrast'; syncControls(); state.notice = 'The teeth model couldn’t run here, so teeth were separated by brightness instead.'; }
  else {
    S.teeth = { prob: m.prob, w: m.w, h: m.h }; S.teethMs = m.ms; buildTeethOverlay();
    if (worker) worker.postMessage({ type: 'teeth', sid: S.sid, prob: m.prob, w: m.w, h: m.h });
  }
  setStage();
  if (S.after) { const f = S.after; S.after = null; f(); }
  else if (m.error && state.notice) { setStatus(state.notice, true); state.notice = null; }
  request('opg');
}
function buildTeethOverlay() {
  const S = state.src, t = S.teeth; if (!t) return;
  const c = S.teethCanvas || (S.teethCanvas = document.createElement('canvas'));
  c.width = t.w; c.height = t.h;
  const cx = c.getContext('2d'), id = cx.createImageData(t.w, t.h), d = id.data;
  for (let i = 0, j = 0; i < t.prob.length; i++, j += 4) { const a = smoothstep(0.35, 0.65, t.prob[i]); d[j] = 229; d[j + 1] = 96; d[j + 2] = 80; d[j + 3] = a * 150; }
  cx.putImageData(id, 0, 0);
}
function requestBuild(reset) {
  pending = { params: { ...state.params }, thr: state.render.thr, reset: !!reset || !!(pending && pending.reset) };
  if (!current) nextJob();
}
function nextJob() {
  const job = pending; pending = null;
  if (!job || !state.src) { showBusy(false); setStage(); return; }
  job.id = ++buildSeq; current = job; showBusy(true); setStage();
  if (!state.vol) setStatus('Building the 3D model…');
  runJob(job);
}
function runJob(job) {
  if (worker) { sendSourceOnce(); worker.postMessage({ type: 'build', id: job.id, params: job.params, thr: job.thr }); return; }
  setTimeout(() => {
    if (current !== job) return;
    try {
      const S = state.src, r = SH.buildVolume(S.norm, S.prep, S.w, S.h, job.params, mainCache, job.params.teethSrc === 'model' ? S.teeth : null);
      r.mesh = SH.buildMeshes(r.rgba, r.nx, r.ny, r.nz, r.bmin, r.vox, job.thr, r.cmin, r.cmax);
      finishBuild(r, job);
    } catch (e) { console.error(e); setStatus('The reconstruction failed on this image. Try another export of the OPG.', true); current = null; showBusy(false); setStage(); }
  }, 16);
}
let sentNorm = null;
function sendSourceOnce() { const S = state.src; if (S && S.norm !== sentNorm) { sendSource(); sentNorm = S.norm; } }
function finishBuild(r, job) {
  current = null;
  const first = !state.vol || job.reset;
  state.vol = r; state.mesh = r.mesh; r.mesh = null;
  state.arch = SH.buildArch(r.archW, r.archD);
  uploadMeshes(); softGeo = null;
  if (first) defaultCrosshair(); else clampP();
  const S = state.src, tsrc = job.params.teethSrc === 'model' && S.teeth ? ` · teeth by U-Net (${(S.teethMs / 1000).toFixed(1)} s)` : '';
  const tris = state.mesh ? state.mesh.hard.fc + state.mesh.teeth.fc + state.mesh.gum.fc : 0;
  setStatus(`${r.nx}×${r.ny}×${r.nz} voxels · ${r.vox.toFixed(2)} mm · ${(tris / 1000).toFixed(0)}k triangles · ${Math.round(r.ms + (state.mesh ? state.mesh.ms : 0))} ms${tsrc}`);
  if (state.notice) { setStatus(state.notice, true); state.notice = null; }
  request('gl', 'panes', 'opg', 'readout');
  if (pending) nextJob();
  else { showBusy(false); setStage(); if (state.mesh && Math.abs(state.mesh.thr - state.render.thr) > 1e-6) requestMesh(); }
}
/* Re-extract the surfaces for a new threshold without rebuilding the volume. */
let meshSeq = 0, meshBusy = false;
function requestMesh() {
  if (!state.vol || current || meshBusy) return;
  const thr = state.render.thr, id = ++meshSeq, V = state.vol;
  if (worker) { meshBusy = true; worker.postMessage({ type: 'mesh', id, thr, vms: V.ms }); return; }
  meshBusy = true;
  setTimeout(() => {
    meshBusy = false;
    if (state.vol !== V) return;
    state.mesh = SH.buildMeshes(V.rgba, V.nx, V.ny, V.nz, V.bmin, V.vox, thr, V.cmin, V.cmax);
    uploadMeshes(); softGeo = null; request('gl');
    if (Math.abs(state.mesh.thr - state.render.thr) > 1e-6) requestMesh();
  }, 16);
}
function onMeshed(m) {
  meshBusy = false;
  if (m.error) { console.error(m.error); return; }
  if (!state.vol || m.vms !== state.vol.ms || current) return;
  state.mesh = m.mesh; uploadMeshes(); softGeo = null; request('gl');
  if (Math.abs(m.mesh.thr - state.render.thr) > 1e-6) requestMesh();
}
function showBusy(on) { $('#busy').hidden = !on; }
function setStatus(msg, warn) { const el = $('#status'); el.textContent = msg; el.title = msg; el.classList.toggle('warn', !!warn); }
/* Small progress chip in the 3D view: first build, then the U-Net refinement running in the background. */
function setStage() {
  const el = $('#gl-stage'); if (!el) return;
  const S = state.src;
  let msg = '';
  if (!state.vol && (current || pending)) msg = 'Building the 3D model…';
  else if (S && S.teethPending && wantModel()) msg = 'Refining the teeth with the trained model…';
  else if (current) msg = 'Updating the 3D model…';
  el.textContent = msg; el.hidden = !msg;
}

/* =====================================================================
   Sources: demo, images, DICOM
   ===================================================================== */
let sidSeq = 0;
function setSource(raw, w, h, name, demo) {
  state.src = { raw, w, h, name, demo, sid: ++sidSeq };
  prepSource(false);
  $('#file-name').textContent = name; $('#file-name').title = name;
  state.measures = [];
}
function prepSource() {
  const S = state.src;
  let g = S.raw;
  if (state.params.invert) { g = new Float32Array(S.raw.length); for (let i = 0; i < g.length; i++) g[i] = 1 - S.raw[i]; }
  S.prep = g;
  S.norm = SH.percentileNormalize(g);
  const c = S.canvas || (S.canvas = document.createElement('canvas'));
  c.width = S.w; c.height = S.h;
  const cx = c.getContext('2d'), id = cx.createImageData(S.w, S.h), d = id.data;
  for (let i = 0, j = 0; i < S.norm.length; i++, j += 4) { const v = S.norm[i] * 255; d[j] = d[j + 1] = d[j + 2] = v; d[j + 3] = 255; }
  cx.putImageData(id, 0, 0);
}
function loadDemo() {
  setStatus('Drawing the demo OPG…'); showBusy(true);
  setTimeout(() => {
    const raw = makeDemoOPG();
    Object.assign(state.params, DEFAULT_PARAMS);
    syncControls();
    setSource(raw, DEMO.W, DEMO.H, 'Demo OPG (synthetic)', true);
    // Show the demo straight away (brightness contrast), then swap in the U-Net's teeth once it has run.
    const S = state.src;
    requestBuild(true);
    if (wantModel()) { S.after = () => { if (state.src === S && state.params.teethSrc === 'model' && S.teeth) requestBuild(false); }; requestSegment(); }
  }, 20);
}
async function loadFile(file) {
  if (!file) return;
  try {
    const head = new Uint8Array(await file.slice(0, 132).arrayBuffer());
    const isDicom = /\.dcm$/i.test(file.name || '') || file.type === 'application/dicom' || (head.length >= 132 && head[128] === 68 && head[129] === 73 && head[130] === 67 && head[131] === 77);
    if (isDicom) { setStatus('Reading DICOM…'); const r = await readDicom(file); ingestGray(r.raw, r.w, r.h, file.name || 'DICOM'); return; }
    if (!file.type || !file.type.startsWith('image/')) { setStatus('That file isn’t an image. Load a PNG, JPEG or DICOM OPG.', true); return; }
    const img = await decodeImage(file);
    const { raw, w, h } = imageToGray(img);
    ingestGray(raw, w, h, file.name || 'Pasted image');
  } catch (e) {
    if (!(e && e.userMessage)) console.error(e);
    setStatus(e && e.userMessage ? e.userMessage : 'Couldn’t read that file. Export the OPG again as PNG or JPEG and retry.', true);
  }
}
function decodeImage(blob) {
  return new Promise((res, rej) => {
    const url = URL.createObjectURL(blob), img = new Image();
    img.onload = () => { URL.revokeObjectURL(url); res(img); };
    img.onerror = () => { URL.revokeObjectURL(url); rej(new Error('decode')); };
    img.src = url;
  });
}
function imageToGray(img) {
  const sc = Math.min(1, (LITE ? 1200 : 1600) / img.naturalWidth, (LITE ? 750 : 1000) / img.naturalHeight);
  const w = Math.max(1, Math.round(img.naturalWidth * sc)), h = Math.max(1, Math.round(img.naturalHeight * sc));
  const c = document.createElement('canvas'); c.width = w; c.height = h;
  const cx = c.getContext('2d', { willReadFrequently: true }); cx.drawImage(img, 0, 0, w, h);
  const d = cx.getImageData(0, 0, w, h).data, raw = new Float32Array(w * h);
  for (let i = 0, j = 0; i < raw.length; i++, j += 4) raw[i] = (0.299 * d[j] + 0.587 * d[j + 1] + 0.114 * d[j + 2]) / 255;
  return { raw, w, h };
}
function downscaleGray(raw, w, h, maxW, maxH) {
  const f = Math.max(1, Math.ceil(Math.max(w / maxW, h / maxH)));
  if (f === 1) return { raw, w, h };
  const W = Math.floor(w / f), H = Math.floor(h / f), out = new Float32Array(W * H), inv = 1 / (f * f);
  for (let y = 0; y < H; y++) for (let x = 0; x < W; x++) {
    let s = 0; for (let yy = 0; yy < f; yy++) { const o = (y * f + yy) * w + x * f; for (let xx = 0; xx < f; xx++) s += raw[o + xx]; }
    out[y * W + x] = s * inv;
  }
  return { raw: out, w: W, h: H };
}
function ingestGray(raw, w, h, name) {
  Object.assign(state.params, DEFAULT_PARAMS, { occl: 0.5, curve: 0.05 });
  setSource(raw, w, h, name, false);
  const S = state.src, aspect = w / h;
  const odd = aspect < 1.45 || w < 500;
  const oddMsg = `This image is ${w}\u00d7${h}. OPGs are usually about twice as wide as tall, so the 3D model may not make sense.`;
  const unsure = 'Auto-fit wasn\u2019t confident on this image. Check that the tooth ticks sit on the teeth and adjust Arch fit if they don\u2019t.';
  // 1) quick fit from brightness and an immediate build, so the 3D shows within a second or two
  const fit = autoFit(true); syncControls();
  if (odd) state.notice = oddMsg;
  requestBuild(true);
  // 2) the trained model finds the teeth in the background; refit from them and rebuild
  if (wantModel()) {
    S.after = () => {
      if (state.src !== S || state.params.teethSrc !== 'model' || !S.teeth) return;
      if (!S.fitTouched) { const f2 = fitArch(true); syncControls(); if (!odd && f2 && !f2.ok) state.notice = unsure; }
      requestBuild(false);
    };
    requestSegment();
  } else if (!odd && fit && !fit.ok) state.notice = unsure;
}
/* Arch fit from the segmented teeth: midline and span from the dentition's extent, occlusal curve from the gap between the rows. */
function autoFitFromTeeth(t) {
  const { prob, w, h } = t, col = new Float32Array(w);
  let total = 0;
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) if (prob[y * w + x] > 0.5) { col[x]++; total++; }
  if (total < 0.01 * w * h) return null;
  const cs = new Float32Array(w); let mx = 0;
  for (let x = 0; x < w; x++) { let s = 0, n = 0; for (let k = -3; k <= 3; k++) { const q = x + k; if (q >= 0 && q < w) { s += col[q]; n++; } } cs[x] = s / n; mx = Math.max(mx, cs[x]); }
  const thr = Math.max(3, 0.1 * mx); let L = -1, R = -1;
  for (let x = 0; x < w; x++) if (cs[x] > thr) { if (L < 0) L = x; R = x; }
  if (R - L < w * 0.2) return null;
  const mid = (L + R + 1) / 2 / w, span = clamp((R - L + 1) / w / 0.71, 0.7, 1.3);
  const pts = [], r0 = Math.round(h * 0.18), r1 = Math.round(h * 0.82), reach = Math.round(h * 0.16);
  for (let x = Math.round(L + 0.08 * (R - L)); x <= R - 0.08 * (R - L); x += 4) {
    const prof = new Float32Array(h);
    for (let y = 0; y < h; y++) { let s = 0; for (let k = 0; k < 4 && x + k < w; k++) s += prob[y * w + x + k]; prof[y] = s / 4; }
    let best = Infinity, by = -1;
    for (let y = r0; y <= r1; y++) {
      let up = 0, dn = 0;
      for (let k = 3; k <= reach; k++) { if (y - k >= 0) up = Math.max(up, prof[y - k]); if (y + k < h) dn = Math.max(dn, prof[y + k]); }
      if (up < 0.5 || dn < 0.5) continue;
      if (prof[y] < best - 1e-4) { best = prof[y]; by = y; }
    }
    if (by < 0 || best > 0.5) continue;
    let a = by, b = by; while (a > r0 && prof[a - 1] <= best + 0.05) a--; while (b < r1 && prof[b + 1] <= best + 0.05) b++;
    pts.push({ n: ((x + 2) / w - mid) / (span / 2), y: ((a + b) / 2 + 0.5) / h });
  }
  if (pts.length < 6) return { ok: false, mid, span };
  let fa = 0.5, fc = 0.05, used = pts;
  for (let round = 0; round < 3; round++) {
    let S0 = 0, S1 = 0, S2 = 0, T0 = 0, T1 = 0;
    for (const p of used) { const q = -p.n * p.n; S0++; S1 += q; S2 += q * q; T0 += p.y; T1 += p.y * q; }
    const det = S0 * S2 - S1 * S1; if (Math.abs(det) > 1e-9) { fa = (T0 * S2 - T1 * S1) / det; fc = (S0 * T1 - S1 * T0) / det; }
    const keep = pts.filter(p => Math.abs(fa - fc * p.n * p.n - p.y) < 0.025);
    if (keep.length < 6) break; used = keep;
  }
  return { ok: fa > 0.3 && fa < 0.72, mid, span, occl: clamp(fa, 0.3, 0.72), curve: clamp(fc, -0.2, 0.4), n: used.length };
}
function fitArch(silent) {
  const S = state.src; if (!S) return null;
  const P = state.params;
  const t = S.teeth && P.teethSrc === 'model' ? autoFitFromTeeth(S.teeth) : null;
  if (t && t.ok) { P.mid = clamp(t.mid, 0.38, 0.62); P.span = t.span; P.occl = t.occl; P.curve = t.curve; if (!silent) syncControls(); return { ok: true, mid: P.mid, occl: P.occl, curve: P.curve, span: P.span, by: 'teeth' }; }
  const r = autoFit(silent);
  if (t && t.mid) { P.mid = clamp(t.mid, 0.38, 0.62); P.span = t.span; if (!silent) syncControls(); }
  return r;
}

/* Minimal DICOM reader: uncompressed (implicit/explicit little endian) and JPEG baseline. */
async function readDicom(file) {
  const buf = await file.arrayBuffer(), dv = new DataView(buf), u8 = new Uint8Array(buf), n = buf.byteLength;
  const fail = msg => { const e = new Error(msg); e.userMessage = msg; throw e; };
  const LONG = new Set(['OB', 'OD', 'OF', 'OL', 'OV', 'OW', 'SQ', 'SV', 'UC', 'UN', 'UR', 'UT', 'UV']);
  const str = (o, len) => { let s = ''; for (let i = 0; i < len; i++) { const c = u8[o + i]; if (c) s += String.fromCharCode(c); } return s.trim(); };
  let off = (n > 132 && str(128, 4) === 'DICM') ? 132 : 0;
  let ts = '1.2.840.10008.1.2.1', explicit = true, inMeta = true;
  if (!off && !/^[A-Z]{2}$/.test(String.fromCharCode(u8[4], u8[5]))) ts = '1.2.840.10008.1.2'; // raw dataset without a header: guess implicit VR
  const tags = new Map(), WANT = new Set([0x00280002, 0x00280004, 0x00280006, 0x00280008, 0x00280010, 0x00280011, 0x00280100, 0x00280101, 0x00280103, 0x00281052, 0x00281053]);
  const el = (o, exp) => {
    const group = dv.getUint16(o, true), elem = dv.getUint16(o + 2, true);
    if (group === 0xFFFE) return { group, elem, vr: '', len: dv.getUint32(o + 4, true), data: o + 8 };
    if (exp) {
      const vr = String.fromCharCode(u8[o + 4], u8[o + 5]);
      if (LONG.has(vr)) return { group, elem, vr, len: dv.getUint32(o + 8, true), data: o + 12 };
      return { group, elem, vr, len: dv.getUint16(o + 6, true), data: o + 8 };
    }
    return { group, elem, vr: '', len: dv.getUint32(o + 4, true), data: o + 8 };
  };
  const skipSeq = (o, exp) => {
    while (o + 8 <= n) {
      const g = dv.getUint16(o, true), e = dv.getUint16(o + 2, true), len = dv.getUint32(o + 4, true);
      if (g === 0xFFFE && e === 0xE0DD) return o + 8;
      if (g !== 0xFFFE || e !== 0xE000) return o;
      o = len === 0xFFFFFFFF ? skipItem(o + 8, exp) : o + 8 + len;
    }
    return o;
  };
  const skipItem = (o, exp) => {
    while (o + 8 <= n) {
      if (dv.getUint16(o, true) === 0xFFFE && dv.getUint16(o + 2, true) === 0xE00D) return o + 8;
      const x = el(o, exp);
      o = x.len === 0xFFFFFFFF ? skipSeq(x.data, exp) : x.data + x.len;
    }
    return o;
  };
  let pix = null;
  while (off + 8 <= n) {
    const g = dv.getUint16(off, true);
    if (inMeta && g !== 0x0002) {
      inMeta = false;
      if (ts === '1.2.840.10008.1.2.2') fail('Big-endian DICOM isn’t supported. Export the OPG as PNG or JPEG.');
      if (ts === '1.2.840.10008.1.2.1.99') fail('Deflated DICOM isn’t supported. Export the OPG as PNG or JPEG.');
      explicit = ts !== '1.2.840.10008.1.2';
    }
    const x = el(off, inMeta ? true : explicit);
    const key = ((x.group << 16) | x.elem) >>> 0;
    if (key === 0x7FE00010) { pix = x; break; }
    if (x.len === 0xFFFFFFFF) { off = skipSeq(x.data, inMeta ? true : explicit); continue; }
    if (key === 0x00020010) ts = str(x.data, x.len).replace(/\0/g, '');
    else if (WANT.has(key)) tags.set(key, x);
    off = x.data + x.len;
  }
  if (!pix) fail('This DICOM file has no image data.');
  const us = k => { const t = tags.get(k); return t ? dv.getUint16(t.data, true) : undefined; };
  const ds = k => { const t = tags.get(k); return t ? parseFloat(str(t.data, t.len).split('\\')[0]) : undefined; };
  const cs = k => { const t = tags.get(k); return t ? str(t.data, t.len) : ''; };
  const rows = us(0x00280010), cols = us(0x00280011);
  const bits = us(0x00280100) || 8, stored = us(0x00280101) || bits, signed = us(0x00280103) === 1;
  const spp = us(0x00280002) || 1, planar = us(0x00280006) || 0, photo = cs(0x00280004) || 'MONOCHROME2';
  const slope = ds(0x00281053) ?? 1, icpt = ds(0x00281052) ?? 0;
  if (!rows || !cols) fail('This DICOM file has no image size.');
  const JPEG = ['1.2.840.10008.1.2.4.50', '1.2.840.10008.1.2.4.51'];
  if (pix.len === 0xFFFFFFFF) {
    if (!JPEG.includes(ts)) {
      const kind = /^1\.2\.840\.10008\.1\.2\.4\.(9\d)$/.test(ts) ? 'JPEG 2000' : /\.4\.(57|70|80|81)$/.test(ts) ? 'lossless JPEG' : ts === '1.2.840.10008.1.2.5' ? 'RLE' : 'this compression';
      fail(`This DICOM uses ${kind}, which the browser can’t decode. Export the OPG as PNG or JPEG.`);
    }
    let o = pix.data; const frags = [];
    while (o + 8 <= n) {
      const g = dv.getUint16(o, true), e = dv.getUint16(o + 2, true), len = dv.getUint32(o + 4, true);
      if (g === 0xFFFE && e === 0xE0DD) break;
      if (g !== 0xFFFE || e !== 0xE000) break;
      frags.push(u8.subarray(o + 8, o + 8 + len)); o += 8 + len;
    }
    const data = frags.slice(1).filter(f => f.length);
    const firstFrame = [];
    for (const f of data) { if (firstFrame.length && f[0] === 0xFF && f[1] === 0xD8) break; firstFrame.push(f); }
    let img;
    try { img = await decodeImage(new Blob(firstFrame, { type: 'image/jpeg' })); }
    catch { fail('The JPEG inside this DICOM couldn’t be decoded (12-bit JPEG isn’t supported by browsers). Export the OPG as PNG or JPEG.'); }
    const r = imageToGray(img);
    if (photo === 'MONOCHROME1') for (let i = 0; i < r.raw.length; i++) r.raw[i] = 1 - r.raw[i];
    return r;
  }
  const count = rows * cols, raw = new Float32Array(count);
  const bpp = bits / 8, o0 = pix.data;
  if (o0 + count * spp * bpp > n) fail('This DICOM file is truncated.');
  const mask = stored < 16 ? (1 << stored) - 1 : 0xFFFF;
  const read = i => {
    if (bits === 8) return u8[o0 + i];
    const v = signed ? dv.getInt16(o0 + i * 2, true) : (dv.getUint16(o0 + i * 2, true) & mask);
    return v;
  };
  if (spp === 3) {
    for (let i = 0; i < count; i++) {
      const r = planar ? read(i) : read(i * 3), g = planar ? read(i + count) : read(i * 3 + 1), b = planar ? read(i + 2 * count) : read(i * 3 + 2);
      raw[i] = 0.299 * r + 0.587 * g + 0.114 * b;
    }
  } else for (let i = 0; i < count; i++) raw[i] = read(i) * slope + icpt;
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < count; i++) { const v = raw[i]; if (v < lo) lo = v; if (v > hi) hi = v; }
  const sc = 1 / ((hi - lo) || 1);
  for (let i = 0; i < count; i++) raw[i] = (raw[i] - lo) * sc;
  if (photo === 'MONOCHROME1') for (let i = 0; i < count; i++) raw[i] = 1 - raw[i];
  return downscaleGray(raw, cols, rows, 1600, 1000);
}

/* Auto-fit: midline by left/right symmetry, occlusal curve from the dark gap between the tooth rows. */
function autoFit(silent) {
  const S = state.src; if (!S) return null;
  const f = Math.max(1, Math.round(S.w / 500)), w = Math.floor(S.w / f), h = Math.floor(S.h / f);
  const a = new Float32Array(w * h);
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    let s = 0; for (let yy = 0; yy < f; yy++) for (let xx = 0; xx < f; xx++) s += S.norm[(y * f + yy) * S.w + x * f + xx];
    a[y * w + x] = s / (f * f);
  }
  const b = SH.gaussBlur(a, w, h, 1.5);
  let best = -Infinity, bc = w / 2;
  const y0 = Math.floor(h * 0.25), y1 = Math.floor(h * 0.85);
  for (let c = Math.round(w * 0.4); c <= Math.round(w * 0.6); c++) {
    const half = Math.min(c, w - 1 - c, Math.round(w * 0.3)); let sum = 0, cnt = 0;
    for (let y = y0; y < y1; y += 2) for (let dx = 4; dx < half; dx += 2) { sum -= Math.abs(b[y * w + c - dx] - b[y * w + c + dx]); cnt++; }
    const sc = sum / cnt; if (sc > best) { best = sc; bc = c; }
  }
  const pts = [];
  const k0 = Math.max(1, Math.round(h * 0.012)), k1 = Math.round(h * 0.09);
  for (let bi = -6; bi <= 6; bi++) {
    const cxb = bc + bi * w * 0.045, xa = Math.max(0, Math.round(cxb - w * 0.02)), xb = Math.min(w - 1, Math.round(cxb + w * 0.02));
    const prof = new Float32Array(h);
    for (let y = 0; y < h; y++) { let s = 0; for (let x = xa; x <= xb; x++) s += b[y * w + x]; prof[y] = s / (xb - xa + 1); }
    let bs = -1, by = 0;
    for (let y = Math.round(h * 0.3); y <= Math.round(h * 0.72); y++) {
      let up = 0, dn = 0;
      for (let k = k0; k <= k1; k++) { if (y - k >= 0) up = Math.max(up, prof[y - k]); if (y + k < h) dn = Math.max(dn, prof[y + k]); }
      const sc = Math.min(up, dn) - prof[y];
      if (sc > bs) { bs = sc; by = y; }
    }
    pts.push({ n: (cxb - bc) / (w / 2), y: (by + 0.5) / h, wt: Math.max(0, bs) });
  }
  // Weighted least squares for y = a - c n^2, two rounds with outlier rejection
  let fa = 0.5, fc = 0.05, used = pts;
  for (let round = 0; round < 2; round++) {
    let S0 = 0, S1 = 0, S2 = 0, T0 = 0, T1 = 0;
    for (const p of used) { const q = -p.n * p.n; S0 += p.wt; S1 += p.wt * q; S2 += p.wt * q * q; T0 += p.wt * p.y; T1 += p.wt * p.y * q; }
    const det = S0 * S2 - S1 * S1;
    if (Math.abs(det) > 1e-9) { fa = (T0 * S2 - T1 * S1) / det; fc = (S0 * T1 - S1 * T0) / det; }
    used = pts.filter(p => Math.abs(fa - fc * p.n * p.n - p.y) < 0.03);
    if (used.length < 5) break;
  }
  const meanScore = pts.reduce((s, p) => s + p.wt, 0) / pts.length;
  const ok = meanScore > 0.05 && used.length >= 6 && fa > 0.33 && fa < 0.7;
  const P = state.params;
  P.mid = clamp((bc + 0.5) / w, 0.38, 0.62);
  if (ok) { P.occl = clamp(fa, 0.33, 0.7); P.curve = clamp(fc, 0, 0.16); }
  if (!silent) syncControls();
  return { ok, mid: P.mid, occl: P.occl, curve: P.curve, score: meanScore };
}

function zOccAtK(k) {
  const A = state.arch, P = state.params;
  return (1 - occFracAt(P, P.mid + P.span * A.s[k] / (2 * A.L))) * state.vol.Hmm;
}
function defaultCrosshair() {
  const A = state.arch, V = state.vol;
  let s, z;
  if (state.src.demo) {
    const t = UPPER[4]; s = t.f * A.L;                       // the implant at 25
    const vf = DEMO.occFracAtS(t.c) - (3 + (t.ch + 6) * DEMO.ppm) / DEMO.H;
    z = (1 - vf) * V.Hmm;
  } else {
    s = LOWER[5].f * A.L;                                     // the 36 region
    z = zOccAtK(idxAtS(A, s)) - LOWER[5].ch - 3;
  }
  const k = idxAtS(A, s);
  state.P = [A.x[k], A.y[k], z]; clampP();
}
function clampP() {
  const V = state.vol, e = V.vox * 0.5;
  for (let i = 0; i < 3; i++) state.P[i] = clamp(state.P[i], V.bmin[i] + e, V.bmax[i] - e);
}

/* =====================================================================
   Slice geometry
   ===================================================================== */
function planeGeom(key) {
  const V = state.vol, P = state.P;
  const C = [(V.bmin[0] + V.bmax[0]) / 2, (V.bmin[1] + V.bmax[1]) / 2, (V.bmin[2] + V.bmax[2]) / 2], S = V.size;
  if (key === 'axial') return { n: [0, 0, -1], U: [1, 0, 0], V: [0, -1, 0], c: [C[0], C[1], P[2]], eu: S[0], ev: S[1], lab: ['R', 'L', 'A', 'P'], scroll: [0, 0, 1] };
  if (key === 'coronal') return { n: [0, 1, 0], U: [1, 0, 0], V: [0, 0, -1], c: [C[0], P[1], C[2]], eu: S[0], ev: S[2], lab: ['R', 'L', 'S', 'I'], scroll: [0, 1, 0] };
  if (key === 'sagittal') return { n: [1, 0, 0], U: [0, -1, 0], V: [0, 0, -1], c: [P[0], C[1], C[2]], eu: S[1], ev: S[2], lab: ['A', 'P', 'S', 'I'], scroll: [1, 0, 0] };
  if (state.mpr.obl === 'arch') {
    const A = state.arch, k = nearestIdx(A, P[0], P[1]);
    const U = [-A.nx[k], -A.ny[k], 0], Vv = [0, 0, -1];
    return { n: V3.cross(U, Vv), U, V: Vv, c: [A.x[k], A.y[k], C[2]], eu: Math.min(46, S[0]), ev: S[2], lab: ['Buc', 'Lin', 'S', 'I'], scroll: [A.tx[k], A.ty[k], 0], k };
  }
  const yw = state.mpr.yaw * Math.PI / 180, pt = state.mpr.pitch * Math.PI / 180;
  const n = [Math.cos(pt) * Math.sin(yw), Math.cos(pt) * Math.cos(yw), Math.sin(pt)];
  const U = V3.norm(V3.cross(n, [0, 0, 1])), Vv = V3.cross(n, U);
  const c = V3.sub(C, V3.mul(n, V3.dot(V3.sub(C, P), n)));
  const e = Math.max(S[0], S[1], S[2]);
  return { n, U, V: Vv, c, eu: e, ev: e, lab: null, scroll: n };
}

/* =====================================================================
   3D renderer: precomputed surface meshes drawn with WebGL 2 or WebGL 1.
   Small GLSL ES 1.00 shaders, so it compiles instantly and runs on any GPU, phones included.
   ===================================================================== */
const glc = $('#gl'), glov = $('#gl-ov'), govx = glov.getContext('2d');
const FORCE_SOFT = /[?&]soft\b/.test(location.search), FORCE_GL1 = /[?&]webgl1\b/.test(location.search);
const GL_OPTS = { alpha: false, antialias: true, depth: true, stencil: false, premultipliedAlpha: true, preserveDrawingBuffer: false, powerPreference: 'default' };
let gl = null, glVer = 0;
if (!FORCE_SOFT) {
  try { if (!FORCE_GL1) { gl = glc.getContext('webgl2', GL_OPTS); if (gl) glVer = 2; } } catch (e) { gl = null; }
  if (!gl) try { gl = glc.getContext('webgl', GL_OPTS) || glc.getContext('experimental-webgl', GL_OPTS); if (gl) glVer = 1; } catch (e) { gl = null; }
}
let glOK = !!gl, uintIdx = false, slowGL = false, glLosses = 0;
let meshProg, lineProg, bgProg, capProg, bgBuf, lineBuf, capBuf, capUV, capTex, capKey = '';
const U = {}, UL = {}, UB = {}, UC = {}, GPU = {};

const VS_MESH = `attribute vec3 aPos; attribute vec3 aNrm; attribute vec4 aAttr;
uniform mat4 uVP;
varying vec3 vPos; varying vec3 vNrm; varying vec4 vAttr;
void main(){ vPos = aPos; vNrm = aNrm; vAttr = aAttr; gl_Position = uVP * vec4(aPos, 1.0); }`;
const FS_MESH = `#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
varying vec3 vPos; varying vec3 vNrm; varying vec4 vAttr;
uniform vec3 uEye, uKey, uFill, uClipN, uLow, uHigh, uBack, uXCol;
uniform float uClipOn, uClipD, uMode, uAlpha, uWrap, uSpec, uGloss, uGrain, uSSS, uCrownSpec;
uniform vec2 uFog;
vec3 aces(vec3 x){ return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0); }
void main(){
  if (uClipOn > 0.5 && dot(uClipN, vPos) > uClipD) discard;
  vec3 n = normalize(vNrm);
  vec3 V = normalize(uEye - vPos);
  float t = vAttr.x, metal = vAttr.y, ao = vAttr.z, grain = vAttr.w;
  if (uMode > 1.5) {
    float e = 1.0 - abs(dot(n, V));
    float g = uAlpha * (0.16 + 0.84 * e * e) * (0.6 + 0.8 * t) + metal * uAlpha * 1.5;
    gl_FragColor = vec4(uXCol * g, 1.0);
    return;
  }
  bool inside = !gl_FrontFacing;
  if (inside) n = -n;
  vec3 alb = mix(uLow * (1.0 + uGrain * (grain - 0.5)), uHigh, t);
  alb = mix(alb, vec3(0.58, 0.6, 0.64), metal);
  float ndl = dot(n, uKey);
  float diff = max((ndl + uWrap) / (1.0 + uWrap), 0.0);
  float fill = max(dot(n, uFill), 0.0);
  float hemi = n.z * 0.5 + 0.5;
  float nv = max(dot(n, V), 0.0);
  vec3 H = normalize(uKey + V);
  float sp = mix(mix(uSpec, uCrownSpec, t), 1.4, metal), glo = mix(mix(uGloss, 80.0, t), 120.0, metal);
  float spec = sp * pow(max(dot(n, H), 0.0), glo) * smoothstep(-0.1, 0.3, ndl);
  vec3 amb = mix(vec3(0.10, 0.095, 0.09), vec3(0.33, 0.36, 0.41), hemi);
  vec3 c = alb * (amb * ao + 1.6 * diff * vec3(1.0, 0.96, 0.9) * mix(1.0, ao, 0.65) + 0.34 * fill * ao * vec3(0.72, 0.8, 1.0));
  c += uSSS * alb * pow(1.0 - nv, 2.0) * ao;
  c += spec * vec3(1.0, 0.97, 0.93) * mix(1.0, ao, 0.5);
  c += pow(1.0 - nv, 4.0) * 0.16 * ao * vec3(0.6, 0.7, 0.85);
  if (inside && uMode < 0.5) c = uBack * (0.5 + 0.5 * ao);
  float fog = clamp((length(uEye - vPos) - uFog.x) / (uFog.y - uFog.x), 0.0, 1.0);
  c *= 1.0 - 0.28 * fog;
  vec3 col = pow(aces(c), vec3(1.0 / 2.2));
  float a = 1.0;
  if (uMode > 0.5) a = uAlpha * (1.0 - t) * (0.28 + 0.72 * pow(1.0 - nv, 1.5));
  gl_FragColor = vec4(col, a);
}`;
const VS_LINE = `attribute vec3 aPos; uniform mat4 uVP; void main(){ gl_Position = uVP * vec4(aPos, 1.0); }`;
const FS_LINE = `precision mediump float; uniform vec4 uCol; void main(){ gl_FragColor = uCol; }`;
const VS_BG = `attribute vec2 aPos; varying vec2 vN; void main(){ vN = aPos; gl_Position = vec4(aPos, 0.999, 1.0); }`;
const FS_BG = `precision mediump float; varying vec2 vN; uniform vec3 uTop, uBot;
void main(){ float r = length(vN * vec2(0.9, 1.0)); vec3 c = mix(uBot, uTop, clamp(vN.y * 0.5 + 0.6, 0.0, 1.0)); gl_FragColor = vec4(mix(c, uBot * 0.8, smoothstep(0.6, 1.5, r)), 1.0); }`;
const VS_CAP = `attribute vec3 aPos; attribute vec2 aUV; uniform mat4 uVP; varying vec2 vUV; void main(){ vUV = aUV; gl_Position = uVP * vec4(aPos, 1.0); }`;
const FS_CAP = `precision mediump float; varying vec2 vUV; uniform sampler2D uTex;
void main(){ vec4 c = texture2D(uTex, vUV); if (c.a < 0.5) discard; gl_FragColor = vec4(c.rgb, 1.0); }`;

function compile(vs, fs, attrs) {
  const mk = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS) && !gl.isContextLost()) throw new Error(gl.getShaderInfoLog(s)); return s; };
  const p = gl.createProgram(); gl.attachShader(p, mk(gl.VERTEX_SHADER, vs)); gl.attachShader(p, mk(gl.FRAGMENT_SHADER, fs));
  attrs.forEach((a, i) => gl.bindAttribLocation(p, i, a));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS) && !gl.isContextLost()) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
function createGLObjects() {
  uintIdx = glVer === 2 || !!gl.getExtension('OES_element_index_uint');
  meshProg = compile(VS_MESH, FS_MESH, ['aPos', 'aNrm', 'aAttr']);
  lineProg = compile(VS_LINE, FS_LINE, ['aPos']);
  bgProg = compile(VS_BG, FS_BG, ['aPos']);
  capProg = compile(VS_CAP, FS_CAP, ['aPos', 'aUV']);
  for (const n of ['uVP', 'uEye', 'uKey', 'uFill', 'uClipN', 'uLow', 'uHigh', 'uBack', 'uXCol', 'uClipOn', 'uClipD', 'uMode', 'uAlpha', 'uWrap', 'uSpec', 'uGloss', 'uGrain', 'uSSS', 'uCrownSpec', 'uFog']) U[n] = gl.getUniformLocation(meshProg, n);
  for (const n of ['uVP', 'uCol']) UL[n] = gl.getUniformLocation(lineProg, n);
  for (const n of ['uTop', 'uBot']) UB[n] = gl.getUniformLocation(bgProg, n);
  for (const n of ['uVP', 'uTex']) UC[n] = gl.getUniformLocation(capProg, n);
  bgBuf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, bgBuf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
  lineBuf = gl.createBuffer(); capBuf = gl.createBuffer(); capUV = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, capUV); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1]), gl.STATIC_DRAW);
  capTex = gl.createTexture(); capKey = '';
  for (const k of Object.keys(GPU)) delete GPU[k];
}
function initGL() {
  if (!gl) { enterSoft(); return; }
  try {
    const ext = gl.getExtension('WEBGL_debug_renderer_info');
    const r = String(ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER));
    slowGL = /SwiftShader|llvmpipe|softpipe|Basic Render|Software/i.test(r) && !/[?&]full\b/.test(location.search);
  } catch (e) { /* renderer string hidden */ }
  try { createGLObjects(); }
  catch (e) { console.error(e); enterSoft('Simplified 3D: the graphics driver refused the 3D renderer.'); return; }
  glc.addEventListener('webglcontextlost', e => {
    e.preventDefault(); glOK = false; glLosses++;
    if (glLosses >= 3) enterSoft('Simplified 3D: the graphics card keeps resetting. Reload the page to try again.');
    else showNote('The graphics card reset. Restoring the 3D view…');
  });
  glc.addEventListener('webglcontextrestored', () => {
    if (soft) return;
    try { createGLObjects(); glOK = true; if (state.mesh) uploadMeshes(); hideNote(); request('gl'); }
    catch (e) { console.error(e); enterSoft(); }
  });
}
function showNote(msg) { const n = $('#gl-note'); n.textContent = msg; n.hidden = false; }
function hideNote() { $('#gl-note').hidden = true; }
function showGLMessage(msg) { const el = $('#gl-msg'); el.textContent = msg; el.hidden = false; }

/* Mesh upload. WebGL 1 without 32-bit indices gets an unindexed copy (rare). */
function expandMesh(m) {
  const n = m.idx.length, pos = new Float32Array(n * 3), nrm = new Int8Array(n * 4), attr = new Uint8Array(n * 4);
  for (let i = 0; i < n; i++) {
    const v = m.idx[i];
    pos[i * 3] = m.pos[v * 3]; pos[i * 3 + 1] = m.pos[v * 3 + 1]; pos[i * 3 + 2] = m.pos[v * 3 + 2];
    for (let c = 0; c < 4; c++) { nrm[i * 4 + c] = m.nrm[v * 4 + c]; attr[i * 4 + c] = m.attr[v * 4 + c]; }
  }
  return { pos, nrm, attr, idx: null };
}
function uploadMeshes() {
  if (!glOK || !state.mesh) return;
  for (const kind of ['hard', 'teeth', 'gum']) {
    const old = GPU[kind];
    if (old) for (const b of [old.pos, old.nrm, old.attr, old.idx]) if (b) gl.deleteBuffer(b);
    let m = state.mesh[kind]; if (!m || !m.idx.length) { GPU[kind] = null; continue; }
    const indexed = uintIdx || m.vc < 65536;
    if (!indexed) m = expandMesh(m);
    const mk = (target, data) => { const b = gl.createBuffer(); gl.bindBuffer(target, b); gl.bufferData(target, data, gl.STATIC_DRAW); return b; };
    GPU[kind] = {
      pos: mk(gl.ARRAY_BUFFER, m.pos), nrm: mk(gl.ARRAY_BUFFER, m.nrm), attr: mk(gl.ARRAY_BUFFER, m.attr),
      idx: indexed ? mk(gl.ELEMENT_ARRAY_BUFFER, uintIdx ? m.idx : new Uint16Array(m.idx)) : null,
      count: indexed ? m.idx.length : m.pos.length / 3, type: uintIdx ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT,
    };
  }
  capKey = '';
}
let attribsOn = 0;
function useAttribs(n) {
  for (let i = 0; i < 3; i++) { if (i < n) gl.enableVertexAttribArray(i); else gl.disableVertexAttribArray(i); }
  attribsOn = n;
}
function drawMesh(kind) {
  const g = GPU[kind]; if (!g) return;
  gl.bindBuffer(gl.ARRAY_BUFFER, g.pos); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, g.nrm); gl.vertexAttribPointer(1, 3, gl.BYTE, true, 4, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, g.attr); gl.vertexAttribPointer(2, 4, gl.UNSIGNED_BYTE, true, 0, 0);
  if (g.idx) { gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, g.idx); gl.drawElements(gl.TRIANGLES, g.count, g.type, 0); }
  else gl.drawArrays(gl.TRIANGLES, 0, g.count);
}

function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
  return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]);
}
function lookAt(e, c, up) {
  const z = V3.norm(V3.sub(e, c)), x = V3.norm(V3.cross(up, z)), y = V3.cross(z, x);
  return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0, -V3.dot(x, e), -V3.dot(y, e), -V3.dot(z, e), 1]);
}
function mul4(a, b) {
  const o = new Float32Array(16);
  for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k]; o[i * 4 + j] = s; }
  return o;
}
function mulVec4(m, v) { const o = [0, 0, 0, 0]; for (let i = 0; i < 4; i++) o[i] = m[i] * v[0] + m[4 + i] * v[1] + m[8 + i] * v[2] + m[12 + i] * v[3]; return o; }
function invert4(a) {
  const [a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33] = a;
  const b00 = a00 * a11 - a01 * a10, b01 = a00 * a12 - a02 * a10, b02 = a00 * a13 - a03 * a10, b03 = a01 * a12 - a02 * a11;
  const b04 = a01 * a13 - a03 * a11, b05 = a02 * a13 - a03 * a12, b06 = a20 * a31 - a21 * a30, b07 = a20 * a32 - a22 * a30;
  const b08 = a20 * a33 - a23 * a30, b09 = a21 * a32 - a22 * a31, b10 = a21 * a33 - a23 * a31, b11 = a22 * a33 - a23 * a32;
  let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
  if (!det) return null; det = 1 / det;
  return new Float32Array([
    (a11 * b11 - a12 * b10 + a13 * b09) * det, (a02 * b10 - a01 * b11 - a03 * b09) * det, (a31 * b05 - a32 * b04 + a33 * b03) * det, (a22 * b04 - a21 * b05 - a23 * b03) * det,
    (a12 * b08 - a10 * b11 - a13 * b07) * det, (a00 * b11 - a02 * b08 + a03 * b07) * det, (a32 * b02 - a30 * b05 - a33 * b01) * det, (a20 * b05 - a22 * b02 + a23 * b01) * det,
    (a10 * b10 - a11 * b08 + a13 * b06) * det, (a01 * b08 - a00 * b10 - a03 * b06) * det, (a30 * b04 - a31 * b02 + a33 * b00) * det, (a21 * b02 - a20 * b04 - a23 * b00) * det,
    (a11 * b07 - a10 * b09 - a12 * b06) * det, (a00 * b09 - a01 * b07 + a02 * b06) * det, (a31 * b01 - a30 * b03 - a32 * b00) * det, (a20 * b03 - a21 * b01 + a22 * b00) * det]);
}
const FOV = 32 * Math.PI / 180;
function camera(aspect) {
  const V = state.vol, c = state.cam;
  const C = V3.add([(V.cmin[0] + V.cmax[0]) / 2, (V.cmin[1] + V.cmax[1]) / 2, (V.cmin[2] + V.cmax[2]) / 2], c.pan);
  const diag = Math.hypot(V.size[0], V.size[1], V.size[2]);
  const rad = Math.hypot(V.cmax[0] - V.cmin[0], V.cmax[1] - V.cmin[1], V.cmax[2] - V.cmin[2]) / 2;
  const fovX = 2 * Math.atan(Math.tan(FOV / 2) * aspect);
  const dist = 0.92 * rad / Math.sin(Math.min(FOV, fovX) / 2) * c.zoom;
  const cp = Math.cos(c.pitch);
  const eye = [C[0] + dist * cp * Math.sin(c.yaw), C[1] + dist * cp * Math.cos(c.yaw), C[2] + dist * Math.sin(c.pitch)];
  const view = lookAt(eye, C, [0, 0, 1]);
  const reach = 2.5 * rad + V3.len(c.pan);
  const vp = mul4(perspective(FOV, aspect, Math.max(dist * 0.02, dist - reach), dist + reach), view);
  const right = [view[0], view[4], view[8]], up = [view[1], view[5], view[9]], back = [view[2], view[6], view[10]];
  const key = V3.norm(V3.add(V3.add(V3.mul(right, -0.5), V3.mul(up, 0.75)), V3.mul(back, 0.55)));
  return { eye, view, vp, inv: invert4(vp), dist, diag, rad, key, C };
}
const CSSV = getComputedStyle(document.documentElement);
const PLANE_COL = { axial: CSSV.getPropertyValue('--ax').trim(), coronal: CSSV.getPropertyValue('--co').trim(), sagittal: CSSV.getPropertyValue('--sa').trim(), oblique: CSSV.getPropertyValue('--ob').trim() };
const hexToRGB = h => { h = (h || '#ffffff').replace('#', ''); return [parseInt(h.slice(0, 2), 16) / 255, parseInt(h.slice(2, 4), 16) / 255, parseInt(h.slice(4, 6), 16) / 255]; };
const PLANE_RGB = Object.fromEntries(Object.entries(PLANE_COL).map(([k, v]) => [k, hexToRGB(v)]));

function clipPlane(M) {
  const r = state.render.clip; if (r === 'none') return null;
  const P = state.P; let n;
  if (r === 'axial') n = [0, 0, 1]; else if (r === 'coronal') n = [0, 1, 0]; else if (r === 'sagittal') n = [1, 0, 0]; else n = planeGeom('oblique').n;
  if (V3.dot(n, V3.sub(M.eye, P)) < 0) n = V3.mul(n, -1);   // always remove the half facing the camera
  return { n, d: V3.dot(n, P) };
}
function sizeGL() {
  const dpr = window.devicePixelRatio || 1;
  const area = Math.max(1, glc.clientWidth * glc.clientHeight);
  const budget = slowGL ? 0.35e6 : LITE ? 1.8e6 : 3e6;
  const sc = Math.min(dpr, 2, Math.sqrt(budget / area));
  const w = Math.max(1, Math.round(glc.clientWidth * sc)), h = Math.max(1, Math.round(glc.clientHeight * sc));
  if (glc.width !== w || glc.height !== h) { glc.width = w; glc.height = h; }
}

/* Materials (linear colour) */
const MAT = {
  bone: [0.58, 0.45, 0.31], enamel: [0.80, 0.77, 0.68], root: [0.70, 0.57, 0.39], gum: [0.50, 0.12, 0.13],
  boneCut: [0.55, 0.44, 0.3], gumCut: [0.42, 0.1, 0.11], xray: [0.8, 0.88, 1.0],
};
const MATP = {
  gum: { low: MAT.gum, high: MAT.gum, back: MAT.gumCut, wrap: 0.55, spec: 0.32, gloss: 55, crown: 0.32, grain: 0.12, sss: 0.55 },
  teeth: { low: MAT.root, high: MAT.enamel, back: MAT.root, wrap: 0.35, spec: 0.16, gloss: 28, crown: 0.85, grain: 0.05, sss: 0.12 },
  hard: { low: MAT.bone, high: MAT.enamel, back: MAT.boneCut, wrap: 0.18, spec: 0.07, gloss: 18, crown: 0.85, grain: 0.32, sss: 0.06 },
};
function setMaterial(kind, mode, alpha) {
  const p = MATP[kind];
  gl.uniform1f(U.uMode, mode); gl.uniform1f(U.uAlpha, alpha);
  gl.uniform3fv(U.uLow, p.low); gl.uniform3fv(U.uHigh, p.high); gl.uniform3fv(U.uBack, p.back);
  gl.uniform1f(U.uWrap, p.wrap); gl.uniform1f(U.uSpec, p.spec); gl.uniform1f(U.uGloss, p.gloss); gl.uniform1f(U.uCrownSpec, p.crown); gl.uniform1f(U.uGrain, p.grain); gl.uniform1f(U.uSSS, p.sss);
}

/* The cut face: the volume sampled on the clip plane, so a cut shows real cross-section structure. */
function capImage(key, gums, boneVis) {
  const V = state.vol, G = planeGeom(key), T = state.render.thr * 255;
  // cover wherever the plane crosses the model: project the model's bounding box onto the plane
  let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
  for (let q = 0; q < 8; q++) {
    const p = [(q & 1 ? V.cmax : V.cmin)[0], (q & 2 ? V.cmax : V.cmin)[1], (q & 4 ? V.cmax : V.cmin)[2]], d = V3.sub(p, G.c);
    const u = V3.dot(d, G.U), v = V3.dot(d, G.V); u0 = Math.min(u0, u); u1 = Math.max(u1, u); v0 = Math.min(v0, v); v1 = Math.max(v1, v);
  }
  u0 -= 2; v0 -= 2; u1 += 2; v1 += 2;
  const res = LITE || slowGL ? 220 : 360, mpp = Math.max(u1 - u0, v1 - v0) / res;
  const w = Math.max(2, Math.ceil((u1 - u0) / mpp)), h = Math.max(2, Math.ceil((v1 - v0) / mpp));
  u1 = u0 + w * mpp; v1 = v0 + h * mpp;
  const img = new Uint8Array(w * h * 4), rg = V.rgba, nx = V.nx, ny = V.ny, nz = V.nz, nxy = nx * ny, iv = 1 / V.vox;
  const smp = new Float32Array(4);
  for (let j = 0; j < h; j++) for (let i = 0; i < w; i++) {
    const du = u0 + (i + 0.5) * mpp, dv = v0 + (j + 0.5) * mpp;
    const gx = (G.c[0] + G.U[0] * du + G.V[0] * dv - V.bmin[0]) * iv - 0.5, gy = (G.c[1] + G.U[1] * du + G.V[1] * dv - V.bmin[1]) * iv - 0.5, gz = (G.c[2] + G.U[2] * du + G.V[2] * dv - V.bmin[2]) * iv - 0.5;
    const o = (j * w + i) * 4;
    if (gx < 0 || gy < 0 || gz < 0 || gx > nx - 1.001 || gy > ny - 1.001 || gz > nz - 1.001) continue;
    const x0 = gx | 0, y0 = gy | 0, z0 = gz | 0, fx = gx - x0, fy = gy - y0, fz = gz - z0, b = (x0 + nx * y0 + nxy * z0) * 4;
    for (let c = 0; c < 4; c++) {
      const q = b + c;
      const c00 = rg[q] + (rg[q + 4] - rg[q]) * fx, c10 = rg[q + nx * 4] + (rg[q + nx * 4 + 4] - rg[q + nx * 4]) * fx;
      const c01 = rg[q + nxy * 4] + (rg[q + nxy * 4 + 4] - rg[q + nxy * 4]) * fx, c11 = rg[q + nxy * 4 + nx * 4] + (rg[q + nxy * 4 + nx * 4 + 4] - rg[q + nxy * 4 + nx * 4]) * fx;
      const c0 = c00 + (c10 - c00) * fy, c1 = c01 + (c11 - c01) * fy; smp[c] = c0 + (c1 - c0) * fz;
    }
    const r = smp[0], tooth = Math.max(smp[1], smp[3]) / 255, metal = smp[3] / 255;
    if (r >= T && (boneVis || tooth >= 0.5)) {
      const d = r / 255;
      let cr, cg, cb;
      if (metal > 0.5) { cr = 190; cg = 195; cb = 204; }
      else if (tooth >= 0.5) { const k = 0.62 + 0.45 * d; cr = 236 * k; cg = 226 * k; cb = 204 * k; }
      else { const k = 0.42 + 0.7 * d; cr = 206 * k; cg = 172 * k; cb = 126 * k; }
      img[o] = Math.min(255, cr); img[o + 1] = Math.min(255, cg); img[o + 2] = Math.min(255, cb); img[o + 3] = 255;
    } else if (gums && smp[2] >= 128) { img[o] = 176; img[o + 1] = 74; img[o + 2] = 78; img[o + 3] = 255; }
  }
  return { img, w, h, G, u0, u1, v0, v1 };
}
function drawCap(gums, boneVis) {
  const R = state.render, key = R.clip, P = state.P, M0 = state.mpr;
  const k = [key, P.join(','), state.vol.ms, R.thr, gums, boneVis, M0.obl, M0.yaw, M0.pitch].join('|');
  if (k !== capKey) {
    const c = capImage(key, gums, boneVis);
    gl.bindTexture(gl.TEXTURE_2D, capTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, c.w, c.h, 0, gl.RGBA, gl.UNSIGNED_BYTE, c.img);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    const G = c.G, at = (u, v) => V3.add(G.c, V3.add(V3.mul(G.U, u), V3.mul(G.V, v)));
    const p00 = at(c.u0, c.v0), p10 = at(c.u1, c.v0), p11 = at(c.u1, c.v1), p01 = at(c.u0, c.v1);
    gl.bindBuffer(gl.ARRAY_BUFFER, capBuf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([...p00, ...p10, ...p11, ...p00, ...p11, ...p01]), gl.DYNAMIC_DRAW);
    capKey = k;
  }
  gl.useProgram(capProg); useAttribs(2);
  gl.bindBuffer(gl.ARRAY_BUFFER, capBuf); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, capUV); gl.vertexAttribPointer(1, 2, gl.FLOAT, false, 0, 0);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, capTex); gl.uniform1i(UC.uTex, 0);
  gl.uniformMatrix4fv(UC.uVP, false, capVP);
  gl.disable(gl.CULL_FACE);
  gl.drawArrays(gl.TRIANGLES, 0, 6);
}
let capVP = null;
function render3D() {
  if (!glOK) { if (soft) renderSoft3D(); return; }
  if (!state.vol || !glc.clientWidth) return;
  if (gl.isContextLost()) return;
  sizeGL();
  const w = glc.width, h = glc.height;
  const M = camera(w / h), R = state.render, xray = R.mode === 'xray';
  gl.viewport(0, 0, w, h);
  gl.depthMask(true); gl.clearColor(0, 0, 0, 1); gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  // background
  gl.disable(gl.DEPTH_TEST); gl.disable(gl.BLEND); gl.disable(gl.CULL_FACE);
  gl.useProgram(bgProg); useAttribs(1);
  gl.bindBuffer(gl.ARRAY_BUFFER, bgBuf); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  if (xray) { gl.uniform3f(UB.uTop, 0.03, 0.035, 0.045); gl.uniform3f(UB.uBot, 0.01, 0.012, 0.016); }
  else { gl.uniform3f(UB.uTop, 0.125, 0.15, 0.17); gl.uniform3f(UB.uBot, 0.035, 0.045, 0.055); }
  gl.drawArrays(gl.TRIANGLES, 0, 3);
  if (state.mesh && GPU.hard !== undefined) {
    const cp = clipPlane(M);
    gl.useProgram(meshProg); useAttribs(3);
    gl.uniformMatrix4fv(U.uVP, false, M.vp);
    gl.uniform3fv(U.uEye, M.eye); gl.uniform3fv(U.uKey, M.key);
    gl.uniform3fv(U.uFill, V3.norm([-M.key[0], -M.key[1], 0.4]));
    gl.uniform3fv(U.uXCol, MAT.xray);
    gl.uniform2f(U.uFog, M.dist - M.rad, M.dist + M.rad);
    gl.uniform1f(U.uClipOn, cp ? 1 : 0); gl.uniform3fv(U.uClipN, cp ? cp.n : [0, 0, 1]); gl.uniform1f(U.uClipD, cp ? cp.d : 0);
    capVP = M.vp;
    if (xray) {
      gl.disable(gl.DEPTH_TEST); gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE);
      setMaterial('hard', 2, 0.2); drawMesh('hard');
      setMaterial('teeth', 2, 0.24); drawMesh('teeth');
      if (R.gums) { setMaterial('gum', 2, 0.035); drawMesh('gum'); }
    } else {
      gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL); gl.disable(gl.BLEND);
      if (cp) gl.disable(gl.CULL_FACE); else { gl.enable(gl.CULL_FACE); gl.cullFace(gl.BACK); }
      const solidBone = R.boneA >= 0.98, gumSolid = R.gumA >= 0.98;
      if (solidBone) { setMaterial('hard', 0, 1); drawMesh('hard'); }
      else { setMaterial('teeth', 0, 1); drawMesh('teeth'); }
      if (cp) { drawCap(R.gums, solidBone || R.boneA >= 0.5); gl.useProgram(meshProg); useAttribs(3); }
      if (R.gums && gumSolid) { if (!cp) gl.enable(gl.CULL_FACE); setMaterial('gum', 0, 1); drawMesh('gum'); }
      // translucent layers: no depth writes, both sides
      gl.disable(gl.CULL_FACE); gl.enable(gl.BLEND); gl.depthMask(false);
      gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      if (!solidBone && R.boneA > 0.01) { setMaterial('hard', 1, R.boneA); drawMesh('hard'); }
      if (R.gums && !gumSolid) { setMaterial('gum', 1, R.gumA); drawMesh('gum'); }
      gl.depthMask(true); gl.disable(gl.BLEND);
    }
  }
  drawLines3D(M);
  drawTriad(M);
}
function drawLines3D(M) {
  const verts = [], groups = [];
  const push = (key, kind, pts) => { groups.push({ key, kind, start: verts.length / 3, count: pts.length }); for (const p of pts) verts.push(p[0], p[1], p[2]); };
  if (state.render.planes) {
    for (const key of PANE_KEYS) {
      const G = planeGeom(key), hu = V3.mul(G.U, G.eu / 2), hv = V3.mul(G.V, G.ev / 2), c = G.c;
      const p = [V3.sub(V3.sub(c, hu), hv), V3.sub(V3.add(c, hu), hv), V3.add(V3.add(c, hu), hv), V3.add(V3.sub(c, hu), hv)];
      push(key, 'fill', [p[0], p[1], p[2], p[0], p[2], p[3]]);
      push(key, 'line', [p[0], p[1], p[1], p[2], p[2], p[3], p[3], p[0]]);
    }
  }
  const P = state.P, r = 5;
  push('sagittal', 'mark', [V3.add(P, [-r, 0, 0]), V3.add(P, [r, 0, 0])]);
  push('coronal', 'mark', [V3.add(P, [0, -r, 0]), V3.add(P, [0, r, 0])]);
  push('axial', 'mark', [V3.add(P, [0, 0, -r]), V3.add(P, [0, 0, r])]);
  gl.useProgram(lineProg); useAttribs(1);
  gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(verts), gl.DYNAMIC_DRAW);
  gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.uniformMatrix4fv(UL.uVP, false, M.vp);
  gl.disable(gl.DEPTH_TEST); gl.disable(gl.CULL_FACE);
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  for (const g of groups) {
    const c = PLANE_RGB[g.key];
    if (g.kind === 'fill') { gl.uniform4f(UL.uCol, c[0], c[1], c[2], 0.035); gl.drawArrays(gl.TRIANGLES, g.start, g.count); }
    else { gl.uniform4f(UL.uCol, c[0], c[1], c[2], g.kind === 'mark' ? 1 : 0.8); gl.drawArrays(gl.LINES, g.start, g.count); }
  }
  gl.disable(gl.BLEND);
}
function drawTriad(M, keep) {
  const w = glov.width, h = glov.height, d = w / Math.max(1, glov.clientWidth);
  if (!keep) govx.clearRect(0, 0, w, h);
  const v = M.view, ox = 34 * d, oy = h - 44 * d, len = 20 * d;
  govx.lineWidth = 1.5 * d; govx.font = `500 ${10 * d}px "IBM Plex Mono", monospace`; govx.textAlign = 'center'; govx.textBaseline = 'middle';
  for (const [ax, lab] of [[[1, 0, 0], 'L'], [[0, 1, 0], 'A'], [[0, 0, 1], 'S']]) {
    const sx = V3.dot(ax, [v[0], v[4], v[8]]), sy = -V3.dot(ax, [v[1], v[5], v[9]]);
    govx.strokeStyle = 'rgba(220,230,234,.55)'; govx.beginPath(); govx.moveTo(ox, oy); govx.lineTo(ox + sx * len, oy + sy * len); govx.stroke();
    govx.fillStyle = 'rgba(220,230,234,.85)'; govx.fillText(lab, ox + sx * (len + 8 * d), oy + sy * (len + 8 * d));
  }
}
/* Click-to-pick: march the same ray on the CPU copy of the volume. */
function pick3D(clientX, clientY) {
  const V = state.vol; if (!V) return null;
  const r = glc.getBoundingClientRect(), M = camera(r.width / r.height);
  const nx = ((clientX - r.left) / r.width) * 2 - 1, ny = 1 - ((clientY - r.top) / r.height) * 2;
  const a = mulVec4(M.inv, [nx, ny, -1, 1]), b = mulVec4(M.inv, [nx, ny, 1, 1]);
  const ro = [a[0] / a[3], a[1] / a[3], a[2] / a[3]], rd = V3.norm(V3.sub([b[0] / b[3], b[1] / b[3], b[2] / b[3]], ro));
  let tn = 0, tf = Infinity;
  for (let i = 0; i < 3; i++) {
    const inv = 1 / rd[i], t0 = (V.bmin[i] - ro[i]) * inv, t1 = (V.bmax[i] - ro[i]) * inv;
    tn = Math.max(tn, Math.min(t0, t1)); tf = Math.min(tf, Math.max(t0, t1));
  }
  const cp = clipPlane(M);
  if (cp) { const dn = V3.dot(cp.n, rd), dr = V3.dot(cp.n, ro) - cp.d; if (Math.abs(dn) > 1e-6) { const tc = -dr / dn; if (dn > 0) tf = Math.min(tf, tc); else tn = Math.max(tn, tc); } }
  if (tf <= tn) return null;
  const thr = state.render.thr * 255, onlyTeeth = state.render.mode !== 'xray' && state.render.boneA < 0.08;
  for (let t = tn; t <= tf; t += V.vox * 0.5) {
    const p = V3.add(ro, V3.mul(rd, t));
    const i = clamp(Math.floor((p[0] - V.bmin[0]) / V.vox), 0, V.nx - 1), j = clamp(Math.floor((p[1] - V.bmin[1]) / V.vox), 0, V.ny - 1), k = clamp(Math.floor((p[2] - V.bmin[2]) / V.vox), 0, V.nz - 1);
    const o = (i + V.nx * (j + V.ny * k)) * 4;
    if (V.rgba[o] >= thr && (!onlyTeeth || V.rgba[o + 1] >= 128 || V.rgba[o + 3] >= 128)) return V3.add(p, V3.mul(rd, V.vox));
  }
  return null;
}

/* =====================================================================
   Software 3D fallback (no WebGL): the same meshes, rasterised on the CPU with a z-buffer
   ===================================================================== */
let soft = false, softCv = null, softCtx = null, softImg = null, softZ = null, softAcc = null, softReason = '', softGeo = null;
const BAYER = [0, 8, 2, 10, 12, 4, 14, 6, 3, 11, 1, 9, 15, 7, 13, 5];
const LUT_G = new Uint8Array(1025); for (let i = 0; i <= 1024; i++) LUT_G[i] = Math.round(255 * Math.pow(i / 1024, 1 / 2.2));
function enterSoft(reason) {
  soft = true; glOK = false;
  if (!softCv) {
    softCv = document.createElement('canvas'); softCv.id = 'gl-soft'; softCv.setAttribute('aria-hidden', 'true');
    glc.parentNode.insertBefore(softCv, glc);             // the WebGL canvas stays on top, invisible, to catch the mouse
    softCv.addEventListener('contextrestored', () => { softImg = null; request('gl'); });
  }
  softCtx = softCv.getContext('2d');
  glc.style.opacity = '0';
  if (!softCtx) { showGLMessage('The 3D view can’t start in this browser. The slice views still work.'); soft = false; return; }
  softReason = reason || 'WebGL is off in this browser, so the 3D is drawn in software. Turn on hardware acceleration in its settings for smoother rotation.';
  showNote(softReason); clearTimeout(enterSoft.t); enterSoft.t = setTimeout(hideNote, 9000);
  request('gl');
}
/* Per-vertex lighting, the same model as the GPU shader. out: r, g, b (0-255) and alpha per vertex. */
function softShade(g, kind, mode, alpha, M, out) {
  const p = MATP[kind], key = M.key, eye = M.eye, fl = V3.norm([-key[0], -key[1], 0.4]), fog0 = M.dist - M.rad, fogR = 2 * M.rad;
  const { pos, nrm, attr, vc } = g;
  const aces = x => { const v = (x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14); return LUT_G[v <= 0 ? 0 : v >= 1 ? 1024 : (v * 1024) | 0]; };
  for (let v = 0; v < vc; v++) {
    const px = pos[v * 3], py = pos[v * 3 + 1], pz = pos[v * 3 + 2];
    let nx = nrm[v * 4] / 127, ny = nrm[v * 4 + 1] / 127, nz = nrm[v * 4 + 2] / 127;
    let vx = eye[0] - px, vy = eye[1] - py, vz = eye[2] - pz; const dist = Math.hypot(vx, vy, vz); vx /= dist; vy /= dist; vz /= dist;
    let nv = nx * vx + ny * vy + nz * vz;
    const t = attr[v * 4] / 255, metal = attr[v * 4 + 1] / 255, ao = attr[v * 4 + 2] / 255, grain = attr[v * 4 + 3] / 255;
    if (mode === 2) { const e = 1 - Math.abs(nv); out[v * 4] = alpha * (0.16 + 0.84 * e * e) * (0.6 + 0.8 * t) + metal * alpha * 1.5; continue; }
    if (mode === 1 && nv < 0) { nx = -nx; ny = -ny; nz = -nz; nv = -nv; }
    nv = nv < 0 ? 0 : nv;
    const gm = 1 + p.grain * (grain - 0.5);
    let ar = p.low[0] * gm + (p.high[0] - p.low[0] * gm) * t, ag = p.low[1] * gm + (p.high[1] - p.low[1] * gm) * t, ab = p.low[2] * gm + (p.high[2] - p.low[2] * gm) * t;
    ar += (0.58 - ar) * metal; ag += (0.6 - ag) * metal; ab += (0.64 - ab) * metal;
    const ndl = nx * key[0] + ny * key[1] + nz * key[2];
    const diff = Math.max((ndl + p.wrap) / (1 + p.wrap), 0), fill = Math.max(nx * fl[0] + ny * fl[1] + nz * fl[2], 0), hemi = nz * 0.5 + 0.5;
    let hx = key[0] + vx, hy = key[1] + vy, hz = key[2] + vz; const hl = Math.hypot(hx, hy, hz) || 1;
    const nh = Math.max((nx * hx + ny * hy + nz * hz) / hl, 0);
    const sp = (p.spec + (p.crown - p.spec) * t) * (1 - metal) + 1.4 * metal, glo = (p.gloss + (80 - p.gloss) * t) * (1 - metal) + 120 * metal;
    const ss = Math.min(1, Math.max(0, (ndl + 0.1) / 0.4)), spec = sp * Math.pow(nh, glo) * ss * ss * (3 - 2 * ss) * (0.5 + 0.5 * ao);
    const df = 1.6 * diff * (1 - 0.65 + 0.65 * ao), fr = (1 - nv) * (1 - nv), rim = fr * fr * 0.16 * ao, sss = p.sss * fr * ao;
    const fog = 1 - 0.28 * Math.min(1, Math.max(0, (dist - fog0) / fogR));
    const ambR = (0.10 + 0.23 * hemi) * ao, ambG = (0.095 + 0.265 * hemi) * ao, ambB = (0.09 + 0.32 * hemi) * ao;
    out[v * 4] = aces((ar * (ambR + df + 0.34 * fill * ao * 0.72 + sss) + spec + rim * 0.6) * fog);
    out[v * 4 + 1] = aces((ag * (ambG + df * 0.96 + 0.34 * fill * ao * 0.8 + sss) + spec * 0.97 + rim * 0.7) * fog);
    out[v * 4 + 2] = aces((ab * (ambB + df * 0.9 + 0.34 * fill * ao + sss) + spec * 0.93 + rim * 0.85) * fog);
    out[v * 4 + 3] = mode === 1 ? alpha * (1 - t) * (0.28 + 0.72 * Math.pow(1 - nv, 1.5)) : 1;
  }
}
/* mode 0 opaque, 1 screen-door translucent, 2 additive (x-ray) */
function softRaster(g, col, S, w, h, D32, Z, mode, cull, clipD, backCol) {
  const idx = g.idx, fc = idx.length / 3, acc = softAcc;
  for (let f = 0; f < fc; f++) {
    const a = idx[f * 3], b = idx[f * 3 + 1], c = idx[f * 3 + 2];
    if (clipD && (clipD[a] > 0 || clipD[b] > 0 || clipD[c] > 0)) continue;
    const x0 = S[a * 3], y0 = S[a * 3 + 1], x1 = S[b * 3], y1 = S[b * 3 + 1], x2 = S[c * 3], y2 = S[c * 3 + 1];
    if (x0 !== x0 || x1 !== x1 || x2 !== x2) continue;             // behind the camera
    const area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    if (area === 0) continue;
    const back = area > 0;
    if (back && cull) continue;
    let minx = Math.max(0, Math.ceil(Math.min(x0, x1, x2) - 0.5)), maxx = Math.min(w - 1, Math.floor(Math.max(x0, x1, x2) - 0.5));
    let miny = Math.max(0, Math.ceil(Math.min(y0, y1, y2) - 0.5)), maxy = Math.min(h - 1, Math.floor(Math.max(y0, y1, y2) - 0.5));
    if (minx > maxx || miny > maxy) continue;
    const z0 = S[a * 3 + 2], z1 = S[b * 3 + 2], z2 = S[c * 3 + 2], inv = 1 / area;
    const flat = back && backCol;
    for (let py = miny; py <= maxy; py++) {
      const qy = py + 0.5;
      for (let px = minx; px <= maxx; px++) {
        const qx = px + 0.5;
        const w0 = ((x1 - qx) * (y2 - qy) - (x2 - qx) * (y1 - qy)) * inv; if (w0 < -1e-7) continue;
        const w1 = ((x2 - qx) * (y0 - qy) - (x0 - qx) * (y2 - qy)) * inv; if (w1 < -1e-7) continue;
        const w2 = 1 - w0 - w1; if (w2 < -1e-7) continue;
        const k = py * w + px;
        if (mode === 2) { acc[k] += w0 * col[a * 4] + w1 * col[b * 4] + w2 * col[c * 4]; continue; }
        const z = w0 * z0 + w1 * z1 + w2 * z2; if (z >= Z[k]) continue;
        if (mode === 1 && (w0 * col[a * 4 + 3] + w1 * col[b * 4 + 3] + w2 * col[c * 4 + 3]) * 16 <= BAYER[((py & 3) << 2) | (px & 3)]) continue;
        Z[k] = z;
        D32[k] = flat ? backCol : (255 << 24) | (((w0 * col[a * 4 + 2] + w1 * col[b * 4 + 2] + w2 * col[c * 4 + 2]) | 0) << 16) | (((w0 * col[a * 4 + 1] + w1 * col[b * 4 + 1] + w2 * col[c * 4 + 1]) | 0) << 8) | ((w0 * col[a * 4] + w1 * col[b * 4] + w2 * col[c * 4]) | 0);
      }
    }
  }
}
function renderSoft3D() {
  const V = state.vol, mesh = state.mesh; if (!V || !mesh || !softCv || !softCv.clientWidth) return;
  const area = softCv.clientWidth * softCv.clientHeight, dpr = window.devicePixelRatio || 1;
  const sc = Math.min(dpr, Math.sqrt((state.moving ? (LITE ? 0.9e5 : 1.4e5) : (LITE ? 2.6e5 : 4.2e5)) / Math.max(1, area)));
  const w = Math.max(1, Math.round(softCv.clientWidth * sc)), h = Math.max(1, Math.round(softCv.clientHeight * sc));
  if (softCv.width !== w || softCv.height !== h) { softCv.width = w; softCv.height = h; }
  if (!softImg || softImg.width !== w || softImg.height !== h) { softImg = softCtx.createImageData(w, h); softZ = new Float32Array(w * h); softAcc = new Float32Array(w * h); }
  const R = state.render, M = camera(w / h), m = M.vp, xray = R.mode === 'xray', cp = clipPlane(M);
  const D32 = new Uint32Array(softImg.data.buffer), Z = softZ;
  for (let y = 0; y < h; y++) {                          // background gradient, matches the GPU renderer
    const t = y / h, r = xray ? 8 : (32 - 23 * t) | 0, g = xray ? 9 : (38 - 27 * t) | 0, b = xray ? 12 : (43 - 29 * t) | 0;
    D32.fill((255 << 24) | (b << 16) | (g << 8) | r, y * w, y * w + w);
  }
  Z.fill(Infinity);
  if (!softGeo || softGeo.mesh !== mesh) softGeo = { mesh };
  const geo = kind => {
    const g = mesh[kind]; if (!g || !g.vc) return null;
    const e = softGeo[kind] || (softGeo[kind] = { S: new Float32Array(g.vc * 3), col: new Float32Array(g.vc * 4), d: new Float32Array(g.vc) });
    const P = g.pos, S = e.S;
    for (let v = 0; v < g.vc; v++) {
      const x = P[v * 3], y = P[v * 3 + 1], z = P[v * 3 + 2];
      const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
      if (cw <= 1e-3) { S[v * 3] = NaN; continue; }
      S[v * 3] = ((m[0] * x + m[4] * y + m[8] * z + m[12]) / cw * 0.5 + 0.5) * w;
      S[v * 3 + 1] = (0.5 - (m[1] * x + m[5] * y + m[9] * z + m[13]) / cw * 0.5) * h;
      S[v * 3 + 2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / cw;
      if (cp) e.d[v] = cp.n[0] * x + cp.n[1] * y + cp.n[2] * z - cp.d;
    }
    return { g, e };
  };
  const pack = c => (255 << 24) | (LUT_G[Math.min(1024, (c[2] * 1024) | 0)] << 16) | (LUT_G[Math.min(1024, (c[1] * 1024) | 0)] << 8) | LUT_G[Math.min(1024, (c[0] * 1024) | 0)];
  const pass = (kind, mode, alpha, cull) => {
    const q = geo(kind); if (!q) return;
    softShade(q.g, kind, mode, alpha, M, q.e.col);
    softRaster(q.g, q.e.col, q.e.S, w, h, D32, Z, mode, cull && !cp, cp ? q.e.d : null, mode === 0 && cp ? pack(V3.mul(MATP[kind].back, 0.85)) : 0);
  };
  if (xray) {
    softAcc.fill(0);
    pass('hard', 2, 0.2, false); pass('teeth', 2, 0.24, false); if (R.gums) pass('gum', 2, 0.035, false);
    const xc = MAT.xray;
    for (let k = 0; k < w * h; k++) {
      const a = softAcc[k]; if (a <= 0) continue;
      const d = D32[k], r = Math.min(255, (d & 255) + a * xc[0] * 255), g = Math.min(255, ((d >> 8) & 255) + a * xc[1] * 255), b = Math.min(255, ((d >> 16) & 255) + a * xc[2] * 255);
      D32[k] = (255 << 24) | (b << 16) | (g << 8) | r;
    }
  } else {
    const solidBone = R.boneA >= 0.98;
    pass(solidBone ? 'hard' : 'teeth', 0, 1, true);
    if (R.gums) { if (R.gumA >= 0.98) pass('gum', 0, 1, true); else pass('gum', 1, R.gumA, false); }
    if (!solidBone && R.boneA > 0.01) pass('hard', 1, R.boneA, false);
  }
  softCtx.putImageData(softImg, 0, 0);
  // slice-plane outlines and the crosshair marker on the overlay
  const d = glov.width / Math.max(1, glov.clientWidth);
  govx.clearRect(0, 0, glov.width, glov.height);
  const proj = p => { const cw = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15]; return cw <= 0 ? null : [((m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12]) / cw * 0.5 + 0.5) * glov.width, (0.5 - (m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13]) / cw * 0.5) * glov.height]; };
  const line = (a, b, col, wd) => { const A = proj(a), B = proj(b); if (!A || !B) return; govx.strokeStyle = col; govx.lineWidth = wd * d; govx.beginPath(); govx.moveTo(A[0], A[1]); govx.lineTo(B[0], B[1]); govx.stroke(); };
  if (R.planes) for (const k of PANE_KEYS) {
    const G = planeGeom(k), hu = V3.mul(G.U, G.eu / 2), hv = V3.mul(G.V, G.ev / 2), c = G.c;
    const q = [V3.sub(V3.sub(c, hu), hv), V3.sub(V3.add(c, hu), hv), V3.add(V3.add(c, hu), hv), V3.add(V3.sub(c, hu), hv)];
    govx.globalAlpha = 0.8; for (let j = 0; j < 4; j++) line(q[j], q[(j + 1) % 4], PLANE_COL[k], 1); govx.globalAlpha = 1;
  }
  const Pc = state.P, r5 = 5;
  line(V3.add(Pc, [-r5, 0, 0]), V3.add(Pc, [r5, 0, 0]), PLANE_COL.sagittal, 1.5);
  line(V3.add(Pc, [0, -r5, 0]), V3.add(Pc, [0, r5, 0]), PLANE_COL.coronal, 1.5);
  line(V3.add(Pc, [0, 0, -r5]), V3.add(Pc, [0, 0, r5]), PLANE_COL.axial, 1.5);
  drawTriad(M, true);
}

/* =====================================================================
   Slice panes (CPU trilinear sampling)
   ===================================================================== */
const views = $('#views');
const panes = PANE_KEYS.map(key => {
  const el = $(`#pane-${key}`), cv = el.querySelector('canvas');
  return { key, el, cv, ctx: cv.getContext('2d'), info: el.querySelector('.vinfo'), color: PLANE_COL[key] };
});
function sampleAt(V, x, y, z) {
  const i = Math.round((x - V.bmin[0]) / V.vox - 0.5), j = Math.round((y - V.bmin[1]) / V.vox - 0.5), k = Math.round((z - V.bmin[2]) / V.vox - 0.5);
  if (i < 0 || j < 0 || k < 0 || i >= V.nx || j >= V.ny || k >= V.nz) return -1;
  return V.slice[i + V.nx * (j + V.ny * k)];
}
function renderPane(p) {
  const V = state.vol, cv = p.cv, ctx = p.ctx, w = cv.width, h = cv.height;
  if (!V || w < 4 || h < 4) return;
  const G = planeGeom(p.key); p.G = G;
  const mpp = Math.max(G.eu / w, G.ev / h); p.mpp = mpp;
  if (!p.img || p.img.width !== w || p.img.height !== h) p.img = ctx.createImageData(w, h);
  const D = p.img.data, vd = V.slice, nx = V.nx, ny = V.ny, nz = V.nz, nxy = nx * ny, iv = 1 / V.vox;
  const lo = state.mpr.lev - state.mpr.win / 2, sc = 255 / state.mpr.win;
  const o0 = [0, 1, 2].map(a => G.c[a] + (0.5 - w / 2) * mpp * G.U[a] + (0.5 - h / 2) * mpp * G.V[a]);
  const gx0 = (o0[0] - V.bmin[0]) * iv - 0.5, gy0 = (o0[1] - V.bmin[1]) * iv - 0.5, gz0 = (o0[2] - V.bmin[2]) * iv - 0.5;
  const ux = G.U[0] * mpp * iv, uy = G.U[1] * mpp * iv, uz = G.U[2] * mpp * iv;
  const vx = G.V[0] * mpp * iv, vy = G.V[1] * mpp * iv, vz = G.V[2] * mpp * iv;
  const xm = nx - 0.5, ym = ny - 0.5, zm = nz - 0.5;
  let o = 0;
  for (let j = 0; j < h; j++) {
    let gx = gx0 + vx * j, gy = gy0 + vy * j, gz = gz0 + vz * j;
    for (let i = 0; i < w; i++, gx += ux, gy += uy, gz += uz, o += 4) {
      if (gx < -0.5 || gy < -0.5 || gz < -0.5 || gx > xm || gy > ym || gz > zm) { D[o] = 10; D[o + 1] = 14; D[o + 2] = 17; D[o + 3] = 255; continue; }
      let x0 = Math.floor(gx), y0 = Math.floor(gy), z0 = Math.floor(gz), fx = gx - x0, fy = gy - y0, fz = gz - z0;
      if (x0 < 0) { x0 = 0; fx = 0; } else if (x0 > nx - 2) { x0 = nx - 2; fx = 1; }
      if (y0 < 0) { y0 = 0; fy = 0; } else if (y0 > ny - 2) { y0 = ny - 2; fy = 1; }
      if (z0 < 0) { z0 = 0; fz = 0; } else if (z0 > nz - 2) { z0 = nz - 2; fz = 1; }
      const b = x0 + nx * y0 + nxy * z0;
      const c00 = vd[b] + (vd[b + 1] - vd[b]) * fx, c10 = vd[b + nx] + (vd[b + nx + 1] - vd[b + nx]) * fx;
      const c01 = vd[b + nxy] + (vd[b + nxy + 1] - vd[b + nxy]) * fx, c11 = vd[b + nxy + nx] + (vd[b + nxy + nx + 1] - vd[b + nxy + nx]) * fx;
      const c0 = c00 + (c10 - c00) * fy, c1 = c01 + (c11 - c01) * fy;
      let gv = ((c0 + (c1 - c0) * fz) - lo) * sc; gv = gv < 0 ? 0 : gv > 255 ? 255 : gv;
      D[o] = D[o + 1] = D[o + 2] = gv; D[o + 3] = 255;
    }
  }
  ctx.putImageData(p.img, 0, 0);
  drawPaneOverlay(p);
}
function toPix(p, pt) { const G = p.G, d = V3.sub(pt, G.c); return [p.cv.width / 2 + V3.dot(d, G.U) / p.mpp, p.cv.height / 2 + V3.dot(d, G.V) / p.mpp]; }
function drawPaneOverlay(p) {
  const ctx = p.ctx, G = p.G, w = p.cv.width, h = p.cv.height, d = w / Math.max(1, p.cv.clientWidth);
  const P = state.P, [cxp, cyp] = toPix(p, P);
  // watermark
  ctx.save(); ctx.translate(w / 2, h / 2); ctx.rotate(-Math.atan2(h, w) * 0.6);
  ctx.font = `600 ${Math.max(12, Math.min(w, h) * 0.075)}px "IBM Plex Mono", monospace`; ctx.fillStyle = 'rgba(232,163,61,.07)'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
  ctx.fillText('ESTIMATED', 0, 0); ctx.restore();
  if (p.key === 'axial') {
    const A = state.arch;
    ctx.save(); ctx.setLineDash([3 * d, 4 * d]); ctx.strokeStyle = 'rgba(237,227,204,.25)'; ctx.lineWidth = d; ctx.beginPath();
    for (let k = 0; k < A.N; k += 6) { const [x, y] = toPix(p, [A.x[k], A.y[k], P[2]]); k ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
    ctx.stroke(); ctx.restore();
  }
  const nA = G.n, dA = V3.dot(nA, G.c), gap = 9 * d, big = 4 * Math.max(w, h);
  ctx.lineWidth = 1.25 * d;
  for (const q of panes) {
    if (q === p) continue;
    const Hq = planeGeom(q.key), nB = Hq.n, dB = V3.dot(nB, Hq.c);
    const dir = V3.cross(nA, nB), dd = V3.dot(dir, dir); if (dd < 1e-6) continue;
    const pt = V3.mul(V3.add(V3.mul(V3.cross(nB, dir), dA), V3.mul(V3.cross(dir, nA), dB)), 1 / dd);
    const [px, py] = toPix(p, pt);
    let du = V3.dot(dir, G.U), dv = V3.dot(dir, G.V); const l = Math.hypot(du, dv); if (l < 1e-6) continue; du /= l; dv /= l;
    const lc = (cxp - px) * du + (cyp - py) * dv;
    ctx.strokeStyle = q.color; ctx.globalAlpha = 0.85; ctx.beginPath();
    ctx.moveTo(px + du * (lc - big), py + dv * (lc - big)); ctx.lineTo(px + du * (lc - gap), py + dv * (lc - gap));
    ctx.moveTo(px + du * (lc + gap), py + dv * (lc + gap)); ctx.lineTo(px + du * (lc + big), py + dv * (lc + big));
    ctx.stroke(); ctx.globalAlpha = 1;
  }
  // measurements lying in this plane
  ctx.font = `500 ${11 * d}px "IBM Plex Mono", monospace`;
  for (const m of state.measures) {
    if (m.key !== p.key || Math.abs(V3.dot(G.n, m.a) - dA) > state.vol.vox * 0.75 || Math.abs(V3.dot(G.n, m.b) - dA) > state.vol.vox * 0.75) continue;
    const [ax, ay] = toPix(p, m.a), [bx, by] = toPix(p, m.b), L = V3.len(V3.sub(m.b, m.a));
    ctx.strokeStyle = '#EDE3CC'; ctx.lineWidth = 1.5 * d; ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by); ctx.stroke();
    ctx.fillStyle = '#EDE3CC'; for (const [x, y] of [[ax, ay], [bx, by]]) { ctx.beginPath(); ctx.arc(x, y, 2.5 * d, 0, 7); ctx.fill(); }
    const label = `≈${L.toFixed(1)} mm`, tw = ctx.measureText(label).width, mx = (ax + bx) / 2, my = (ay + by) / 2 - 12 * d;
    ctx.fillStyle = 'rgba(10,14,17,.85)'; ctx.fillRect(mx - tw / 2 - 4 * d, my - 8 * d, tw + 8 * d, 16 * d);
    ctx.fillStyle = '#EDE3CC'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(label, mx, my);
  }
  ctx.font = `500 ${10.5 * d}px "IBM Plex Mono", monospace`; ctx.fillStyle = 'rgba(220,230,234,.7)'; ctx.textBaseline = 'middle';
  if (G.lab) {
    ctx.textAlign = 'left'; ctx.fillText(G.lab[0], 6 * d, h / 2);
    ctx.textAlign = 'right'; ctx.fillText(G.lab[1], w - 6 * d, h / 2);
    ctx.textAlign = 'center'; ctx.fillText(G.lab[2], w / 2, 30 * d); ctx.fillText(G.lab[3], w / 2, h - 34 * d);
  }
  const bl = 10 / p.mpp, bx = w - 10 * d - bl, by = h - 22 * d;
  ctx.strokeStyle = 'rgba(220,230,234,.7)'; ctx.lineWidth = d; ctx.beginPath();
  ctx.moveTo(bx, by - 3 * d); ctx.lineTo(bx, by); ctx.lineTo(bx + bl, by); ctx.lineTo(bx + bl, by - 3 * d); ctx.stroke();
  ctx.textAlign = 'center'; ctx.fillText('≈10 mm', bx + bl / 2, by + 9 * d);
  let info;
  if (p.key === 'axial') info = `z ${mm(P[2])} mm`;
  else if (p.key === 'coronal') info = `y ${mm(P[1])} mm`;
  else if (p.key === 'sagittal') info = `x ${mm(P[0])} mm`;
  else if (state.mpr.obl === 'arch') { const sv = state.arch.s[G.k]; info = `s ${sv >= 0 ? '+' : '−'}${Math.abs(sv).toFixed(1)} mm`; }
  else info = `${state.mpr.yaw}° / ${state.mpr.pitch}°`;
  p.info.textContent = info;
}

/* =====================================================================
   Source OPG thumbnail with arch-fit overlay
   ===================================================================== */
const opgCv = $('#opg'), opgCtx = opgCv.getContext('2d'), opgTip = $('#opg-tip');
function renderOPG() {
  const S = state.src; if (!S || !state.vol) return;
  const cw = opgCv.clientWidth || 280, ch = cw * S.h / S.w, d = Math.min(window.devicePixelRatio || 1, 2);
  if (opgCv.width !== Math.round(cw * d) || opgCv.height !== Math.round(ch * d)) { opgCv.width = Math.round(cw * d); opgCv.height = Math.round(ch * d); opgCv.style.height = ch + 'px'; }
  const c = opgCtx; c.setTransform(1, 0, 0, 1, 0, 0); c.imageSmoothingQuality = 'high';
  c.drawImage(S.canvas, 0, 0, opgCv.width, opgCv.height);
  if (state.ui.showTeeth && S.teethCanvas && state.params.teethSrc === 'model') c.drawImage(S.teethCanvas, 0, 0, opgCv.width, opgCv.height);
  c.setTransform(d, 0, 0, d, 0, 0);
  const P = state.params, A = state.arch, V = state.vol;
  const occY = u => occFracAt(P, u) * ch;
  c.setLineDash([3, 3]); c.strokeStyle = 'rgba(83,183,214,.7)'; c.lineWidth = 1;
  for (const sg of [-1, 1]) { const x = (P.mid + sg * P.span / 2) * cw; if (x > 0 && x < cw) { c.beginPath(); c.moveTo(x, 0); c.lineTo(x, ch); c.stroke(); } }
  c.strokeStyle = 'rgba(237,227,204,.55)'; c.beginPath();
  for (let i = 0; i <= 60; i++) { const u = i / 60; i ? c.lineTo(u * cw, occY(u)) : c.moveTo(0, occY(0)); }
  c.stroke(); c.setLineDash([]);
  c.font = '500 8.5px "IBM Plex Mono", monospace'; c.textAlign = 'center'; c.lineWidth = 1;
  for (const up of [true, false]) for (const sg of [-1, 1]) {
    const T = up ? UPPER : LOWER;
    for (let n = 0; n < 8; n++) {
      const u = P.mid + P.span * sg * T[n].f / 2, x = u * cw, yo = occY(u);
      c.strokeStyle = 'rgba(237,227,204,.85)'; c.beginPath();
      c.moveTo(x, up ? yo - 2 : yo + 2); c.lineTo(x, up ? yo - 7 : yo + 7); c.stroke();
      if (n === 2 || n === 5) {
        const q = up ? (sg < 0 ? 1 : 2) : (sg < 0 ? 4 : 3);
        c.fillStyle = 'rgba(0,0,0,.6)'; c.fillRect(x - 7, up ? yo - 19 : yo + 9, 14, 10);
        c.fillStyle = '#EDE3CC'; c.fillText(`${q}${n + 1}`, x, up ? yo - 11 : yo + 17);
      }
    }
  }
  c.font = '600 10px "IBM Plex Mono", monospace'; c.fillStyle = 'rgba(237,227,204,.8)';
  c.textAlign = 'left'; c.fillText('R', 5, 12); c.textAlign = 'right'; c.fillText('L', cw - 5, 12);
  const k = nearestIdx(A, state.P[0], state.P[1]);
  const x = (P.mid + P.span * A.s[k] / (2 * A.L)) * cw, y = (1 - state.P[2] / V.Hmm) * ch;
  c.strokeStyle = '#53B7D6'; c.lineWidth = 1.5; c.beginPath(); c.arc(x, y, 5.5, 0, 7); c.stroke();
  c.beginPath(); c.moveTo(x - 10, y); c.lineTo(x - 7, y); c.moveTo(x + 7, y); c.lineTo(x + 10, y); c.moveTo(x, y - 10); c.lineTo(x, y - 7); c.moveTo(x, y + 7); c.lineTo(x, y + 10); c.stroke();
}
function opgPoint(e) {
  const r = opgCv.getBoundingClientRect(), P = state.params, A = state.arch, V = state.vol;
  const u = clamp((e.clientX - r.left) / r.width, 0, 1), v = clamp((e.clientY - r.top) / r.height, 0, 1);
  const s = clamp((u - P.mid) / P.span * 2 * A.L, -A.L, A.L), k = idxAtS(A, s);
  return { pt: [A.x[k], A.y[k], (1 - v) * V.Hmm], x: e.clientX - r.left, y: e.clientY - r.top };
}

/* =====================================================================
   Readout and tooth chart
   ===================================================================== */
function describe(P) {
  const A = state.arch, prm = state.params;
  const k = nearestIdx(A, P[0], P[1]), s = A.s[k];
  const dist = (P[0] - A.x[k]) * A.nx[k] + (P[1] - A.y[k]) * A.ny[k];
  const T = thickAt(s / A.L) * prm.thick;
  const up = P[2] > zOccAtK(k), f = Math.abs(s) / A.L, list = up ? UPPER : LOWER;
  if (k <= 1 || k >= A.N - 2 || Math.abs(dist) > T / 2 + 14) return { fdi: '–', name: 'Outside the modelled jaw', dist };
  if (f > list[7].fEnd + 0.03) return { fdi: '–', name: up ? 'Maxillary tuberosity region' : (f > 0.88 ? 'Condyle and upper ramus' : 'Mandibular ramus'), dist };
  let n = 0, bd = 1e9; list.forEach((t, i) => { const dd = Math.abs(t.f - f); if (dd < bd) { bd = dd; n = i; } });
  const q = up ? (s < 0 ? 1 : 2) : (s < 0 ? 4 : 3);
  return { fdi: `${q}${n + 1}`, name: `${up ? 'Upper' : 'Lower'} ${s < 0 ? 'right' : 'left'} ${TOOTH_NAMES[n]} region`, dist };
}
const chartBtns = new Map();
function buildChart() {
  const chart = $('#chart');
  for (const up of [true, false]) {
    const row = document.createElement('div'); row.className = 'chart-row';
    const left = up ? [18, 17, 16, 15, 14, 13, 12, 11] : [48, 47, 46, 45, 44, 43, 42, 41];
    const right = up ? [21, 22, 23, 24, 25, 26, 27, 28] : [31, 32, 33, 34, 35, 36, 37, 38];
    const mk = fdi => { const b = document.createElement('button'); b.type = 'button'; b.textContent = fdi; b.title = `Jump to ${fdi}`; b.addEventListener('click', () => jumpToTooth(fdi)); chartBtns.set(String(fdi), b); return b; };
    left.forEach(f => row.appendChild(mk(f)));
    row.appendChild(document.createElement('span'));
    right.forEach(f => row.appendChild(mk(f)));
    chart.appendChild(row);
  }
}
function jumpToTooth(fdi) {
  if (!state.vol) return;
  const q = Math.floor(fdi / 10), n = fdi % 10, up = q === 1 || q === 2, sg = (q === 1 || q === 4) ? -1 : 1;
  const T = (up ? UPPER : LOWER)[n - 1], A = state.arch;
  const k = idxAtS(A, sg * T.f * A.L), zo = zOccAtK(k);
  state.P = [A.x[k], A.y[k], zo + (up ? 1 : -1) * (T.ch * 0.9 + 2)]; clampP();
  request('gl', 'panes', 'opg', 'readout');
}
function renderReadout() {
  const r = describe(state.P), P = state.P, v = sampleAt(state.vol, P[0], P[1], P[2]);
  $('#ro-fdi').textContent = r.fdi; $('#ro-name').textContent = r.name; $('#ro-name').title = r.name;
  $('#ro-xyz').innerHTML = `x ${mm(P[0])} · y ${mm(P[1])} · z ${mm(P[2])} mm (est.)<br>${r.dist >= 0 ? 'buccal' : 'lingual'} ${Math.abs(r.dist).toFixed(1)} mm from arch · value ${v < 0 ? '–' : v}`;
  chartBtns.forEach((b, k) => b.classList.toggle('on', k === r.fdi));
}

/* =====================================================================
   Scheduling and sizing
   ===================================================================== */
const dirty = new Set(); let raf = 0;
function request(...what) { what.forEach(w => dirty.add(w)); if (!raf) raf = requestAnimationFrame(frame); }
function frame() {
  raf = 0; if (!state.vol) { dirty.clear(); return; }
  const d = new Set(dirty); dirty.clear();
  if (d.has('gl')) render3D();
  if (d.has('panes')) panes.forEach(renderPane);
  if (d.has('opg')) renderOPG();
  if (d.has('readout')) renderReadout();
}
let settleTimer = 0;
function moving() { state.moving = true; clearTimeout(settleTimer); settleTimer = setTimeout(() => { state.moving = false; request('gl'); }, 180); }
let resizeRaf = 0;
function resizeAll() {
  resizeRaf = 0;
  const dO = Math.min(window.devicePixelRatio || 1, 2);
  glov.width = Math.max(1, Math.round(glov.clientWidth * dO)); glov.height = Math.max(1, Math.round(glov.clientHeight * dO));
  for (const p of panes) {
    const cw = p.cv.clientWidth, chh = p.cv.clientHeight;
    const dp = Math.min(window.devicePixelRatio || 1, Math.sqrt((LITE ? 150000 : 260000) / Math.max(1, cw * chh)), 2);
    p.cv.width = Math.max(1, Math.round(cw * dp)); p.cv.height = Math.max(1, Math.round(chh * dp));
  }
  request('gl', 'panes', 'opg');
}
const scheduleResize = () => { if (!resizeRaf) resizeRaf = requestAnimationFrame(resizeAll); };

/* =====================================================================
   Interaction
   ===================================================================== */
const ptrs = new Map(); let pinch0 = 0, down = null;
glc.addEventListener('pointerdown', e => {
  glc.setPointerCapture(e.pointerId); ptrs.set(e.pointerId, { x: e.clientX, y: e.clientY, pan: e.shiftKey || e.button === 2 });
  if (ptrs.size === 1) down = { x: e.clientX, y: e.clientY, t: performance.now(), moved: 0 };
  if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; pinch0 = Math.hypot(a.x - b.x, a.y - b.y); down = null; }
});
glc.addEventListener('pointermove', e => {
  const p = ptrs.get(e.pointerId); if (!p || !state.vol) return;
  const dx = e.clientX - p.x, dy = e.clientY - p.y; p.x = e.clientX; p.y = e.clientY;
  if (down) down.moved += Math.abs(dx) + Math.abs(dy);
  if (down && down.moved < 4) return;
  const c = state.cam;
  if (ptrs.size === 2) {
    const [a, b] = [...ptrs.values()], dd = Math.hypot(a.x - b.x, a.y - b.y);
    if (pinch0 > 0 && dd > 0) c.zoom = clamp(c.zoom * pinch0 / dd, 0.25, 3); pinch0 = dd;
  } else if (p.pan) {
    const M = camera(glc.clientWidth / glc.clientHeight), v = M.view;
    const k = 2 * M.dist * Math.tan(FOV / 2) / Math.max(1, glc.clientHeight);
    c.pan = V3.add(c.pan, V3.add(V3.mul([v[0], v[4], v[8]], -dx * k), V3.mul([v[1], v[5], v[9]], dy * k)));
  } else { c.yaw -= dx * 0.009; c.pitch = clamp(c.pitch + dy * 0.009, -1.45, 1.45); }
  moving(); request('gl');
});
const endPtr = e => {
  const wasTap = down && ptrs.size === 1 && down.moved < 4 && performance.now() - down.t < 500 && e.type === 'pointerup' && e.button !== 2;
  ptrs.delete(e.pointerId); if (ptrs.size < 2) pinch0 = 0;
  if (wasTap) {
    const hit = pick3D(e.clientX, e.clientY);
    if (hit) { state.P = hit; clampP(); request('gl', 'panes', 'opg', 'readout'); }
  }
  if (!ptrs.size) down = null;
};
glc.addEventListener('pointerup', endPtr); glc.addEventListener('pointercancel', endPtr);
glc.addEventListener('contextmenu', e => e.preventDefault());
glc.addEventListener('wheel', e => { e.preventDefault(); state.cam.zoom = clamp(state.cam.zoom * Math.exp(e.deltaY * 0.0012), 0.25, 3); moving(); request('gl'); }, { passive: false });
glc.addEventListener('keydown', e => {
  const c = state.cam, k = { ArrowLeft: [0.12, 0], ArrowRight: [-0.12, 0], ArrowUp: [0, -0.1], ArrowDown: [0, 0.1] }[e.key];
  if (!k) return; e.preventDefault(); c.yaw += k[0]; c.pitch = clamp(c.pitch + k[1], -1.45, 1.45); moving(); request('gl');
});
function resetView() { Object.assign(state.cam, { yaw: -0.5, pitch: 0.22, zoom: 1, pan: [0, 0, 0] }); request('gl'); }
let spinRaf = 0, spinLast = 0;
function setSpin(on) {
  state.render.spin = on; $('#btn-spin').setAttribute('aria-pressed', String(on));
  if (on && !spinRaf) { spinLast = 0; spinRaf = requestAnimationFrame(spinTick); }
}
function spinTick(ts) {
  if (!state.render.spin) { spinRaf = 0; return; }
  const dt = spinLast ? Math.min(0.05, (ts - spinLast) / 1000) : 0; spinLast = ts;
  state.cam.yaw += dt * 0.45; moving(); request('gl');
  spinRaf = requestAnimationFrame(spinTick);
}

function paneWorld(p, e) {
  const r = p.cv.getBoundingClientRect(), G = p.G; if (!G) return null;
  const px = (e.clientX - r.left) * p.cv.width / r.width, py = (e.clientY - r.top) * p.cv.height / r.height;
  return V3.add(G.c, V3.add(V3.mul(G.U, (px - p.cv.width / 2) * p.mpp), V3.mul(G.V, (py - p.cv.height / 2) * p.mpp)));
}
function stepSlice(p, dir) {
  const G = p.G || planeGeom(p.key);
  state.P = V3.add(state.P, V3.mul(G.scroll, dir * state.vol.vox)); clampP();
  moving(); request('gl', 'panes', 'opg', 'readout');
}
for (const p of panes) {
  p.cv.addEventListener('pointerdown', e => {
    if (!state.vol) return;
    p.cv.setPointerCapture(e.pointerId);
    const w = paneWorld(p, e); if (!w) return;
    if (state.measure) { p.mdrag = { key: p.key, a: w, b: w }; state.measures.push(p.mdrag); request('panes'); return; }
    p.drag = true; state.P = w; if (p.key === 'axial') state.P[2] = p.G.c[2]; clampP(); moving(); request('gl', 'panes', 'opg', 'readout');
  });
  p.cv.addEventListener('pointermove', e => {
    if (p.mdrag) { const w = paneWorld(p, e); if (w) { p.mdrag.b = w; request('panes'); } return; }
    if (p.drag) { const w = paneWorld(p, e); if (w) { state.P = w; if (p.key === 'axial') state.P[2] = p.G.c[2]; clampP(); moving(); request('gl', 'panes', 'opg', 'readout'); } }
  });
  const end = () => { if (p.mdrag && V3.len(V3.sub(p.mdrag.b, p.mdrag.a)) < 0.3) state.measures.pop(); p.mdrag = null; p.drag = false; request('panes'); };
  p.cv.addEventListener('pointerup', end); p.cv.addEventListener('pointercancel', end);
  p.cv.addEventListener('wheel', e => { e.preventDefault(); if (state.vol) stepSlice(p, (e.deltaY < 0 ? 1 : -1) * (e.shiftKey ? 5 : 1)); }, { passive: false });
  p.cv.addEventListener('keydown', e => {
    const m = { ArrowUp: 1, ArrowRight: 1, ArrowDown: -1, ArrowLeft: -1, PageUp: 5, PageDown: -5 }[e.key];
    if (m && state.vol) { e.preventDefault(); stepSlice(p, m); }
  });
}
let opgDrag = false;
opgCv.addEventListener('pointerdown', e => { if (!state.vol) return; opgCv.setPointerCapture(e.pointerId); opgDrag = true; state.P = opgPoint(e).pt; clampP(); moving(); request('gl', 'panes', 'opg', 'readout'); });
opgCv.addEventListener('pointermove', e => {
  if (!state.vol) return;
  const r = opgPoint(e);
  if (opgDrag) { state.P = r.pt; clampP(); moving(); request('gl', 'panes', 'opg', 'readout'); }
  const dsc = describe(r.pt); opgTip.hidden = false; opgTip.textContent = dsc.fdi !== '–' ? `${dsc.fdi} · ${dsc.name.replace(' region', '')}` : dsc.name;
  opgTip.style.left = r.x + 'px'; opgTip.style.top = r.y + 'px';
});
opgCv.addEventListener('pointerleave', () => { opgTip.hidden = true; });
opgCv.addEventListener('pointerup', () => { opgDrag = false; });
opgCv.addEventListener('pointercancel', () => { opgDrag = false; });

// loading
$('#btn-demo').addEventListener('click', loadDemo);
$('#btn-upload').addEventListener('click', () => $('#file').click());
$('#file').addEventListener('change', e => { loadFile(e.target.files[0]); e.target.value = ''; });
const veil = $('#dropveil'); let dragDepth = 0;
const hasFiles = e => e.dataTransfer && [...(e.dataTransfer.types || [])].includes('Files');
window.addEventListener('dragenter', e => { if (!hasFiles(e)) return; e.preventDefault(); dragDepth++; veil.hidden = false; });
window.addEventListener('dragover', e => { if (hasFiles(e)) e.preventDefault(); });
window.addEventListener('dragleave', e => { if (!hasFiles(e)) return; dragDepth = Math.max(0, dragDepth - 1); if (!dragDepth) veil.hidden = true; });
window.addEventListener('drop', e => { e.preventDefault(); dragDepth = 0; veil.hidden = true; if (e.dataTransfer && e.dataTransfer.files[0]) loadFile(e.dataTransfer.files[0]); });
document.addEventListener('paste', e => {
  const item = [...((e.clipboardData && e.clipboardData.items) || [])].find(i => i.type.startsWith('image/'));
  if (item) { e.preventDefault(); loadFile(item.getAsFile()); }
});
$('#btn-autofit').addEventListener('click', () => {
  if (!state.src) return;
  const r = fitArch(false);
  if (r) { state.notice = r.ok ? `Auto-fit${r.by === 'teeth' ? ' from detected teeth' : ''}: midline ${(r.mid * 100).toFixed(1)} %, occlusal plane ${(r.occl * 100).toFixed(1)} % down` : 'Auto-fit wasn’t confident. Set the midline and occlusal plane by hand.'; requestBuild(false); }
});

// presets, layout, maximise
function applyPreset(name) {
  const p = PRESETS[name]; if (!p) return;
  Object.assign(state.render, p, { preset: name });
  syncControls(); syncPresetSeg(); request('gl');
}
function syncPresetSeg() {
  const R = state.render;
  let match = null;
  for (const [k, p] of Object.entries(PRESETS)) if (Object.entries(p).every(([kk, vv]) => R[kk] === vv)) { match = k; break; }
  R.preset = match;
  document.querySelectorAll('#preset-seg button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.preset === match)));
}
document.querySelectorAll('#preset-seg button').forEach(b => b.addEventListener('click', () => applyPreset(b.dataset.preset)));
function setLayout(l) {
  views.dataset.layout = l; setMax('');
  document.querySelectorAll('#layout-seg button').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.layout === l)));
}
document.querySelectorAll('#layout-seg button').forEach(b => b.addEventListener('click', () => setLayout(b.dataset.layout)));
function setMax(key) {
  views.classList.toggle('maxed', !!key);
  document.querySelectorAll('.view').forEach(v => v.classList.toggle('is-max', v.dataset.key === key));
  document.querySelectorAll('.max-btn').forEach(b => { const on = b.dataset.max === key; b.setAttribute('aria-pressed', String(on)); b.title = on ? 'Restore' : 'Maximise'; });
  scheduleResize();
}
document.querySelectorAll('.max-btn').forEach(b => b.addEventListener('click', () => setMax(b.getAttribute('aria-pressed') === 'true' ? '' : b.dataset.max)));
$('#btn-spin').addEventListener('click', () => setSpin(!state.render.spin));
$('#btn-reset').addEventListener('click', resetView);
function setMeasure(on) {
  state.measure = on; views.classList.toggle('measuring', on);
  const b = $('#btn-measure'); if (b) b.setAttribute('aria-pressed', String(on));
}
document.addEventListener('keydown', e => {
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  const t = e.target; if (t && (t.tagName === 'INPUT' || t.tagName === 'SELECT' || t.tagName === 'TEXTAREA')) return;
  const k = e.key.toLowerCase();
  if (k === '1') applyPreset('anatomy'); else if (k === '2') applyPreset('bone'); else if (k === '3') applyPreset('roots'); else if (k === '4') applyPreset('xray');
  else if (k === 'g') { state.render.gums = !state.render.gums; syncControls(); syncPresetSeg(); request('gl'); }
  else if (k === 's') setSpin(!state.render.spin);
  else if (k === 'r') resetView();
  else if (k === 'm') setMeasure(!state.measure);
  else if (e.key === 'Escape') { setMeasure(false); setMax(''); }
});

/* =====================================================================
   Controls
   ===================================================================== */
const controls = [];
function addRange(parent, id, label, obj, key, min, max, step, fmt, onChange) {
  const wrap = document.createElement('div'); wrap.className = 'ctl';
  wrap.innerHTML = `<label for="${id}">${label}</label><output id="${id}-o" for="${id}"></output><input type="range" id="${id}" min="${min}" max="${max}" step="${step}">`;
  parent.appendChild(wrap);
  const inp = wrap.querySelector('input'), out = wrap.querySelector('output');
  const sync = () => { inp.value = obj[key]; out.textContent = fmt(obj[key]); };
  inp.addEventListener('input', () => { obj[key] = parseFloat(inp.value); out.textContent = fmt(obj[key]); onChange(); });
  controls.push(sync); sync();
  return wrap;
}
function addSelect(parent, id, label, obj, key, options, onChange) {
  const wrap = document.createElement('div'); wrap.className = 'ctl';
  wrap.innerHTML = `<label for="${id}">${label}</label><span></span><select id="${id}">${options.map(([v, t]) => `<option value="${v}">${t}</option>`).join('')}</select>`;
  parent.appendChild(wrap);
  const sel = wrap.querySelector('select');
  const sync = () => { sel.value = obj[key]; };
  sel.addEventListener('change', () => { obj[key] = sel.value; onChange(); });
  controls.push(sync); sync();
  return wrap;
}
function addCheck(parent, id, label, obj, key, onChange) {
  const wrap = document.createElement('label'); wrap.className = 'check';
  wrap.innerHTML = `<input type="checkbox" id="${id}"><span>${label}</span>`;
  parent.appendChild(wrap);
  const cb = wrap.querySelector('input');
  const sync = () => { cb.checked = !!obj[key]; };
  cb.addEventListener('change', () => { obj[key] = cb.checked; onChange(); });
  controls.push(sync); sync();
  return wrap;
}
function addRow(parent, html) { const row = document.createElement('div'); row.className = 'row split'; row.innerHTML = html; parent.appendChild(row); return row; }
function syncControls() { controls.forEach(f => f()); }
const pct = v => (v * 100).toFixed(1) + ' %', times = v => '×' + v.toFixed(2);
const P0 = state.params, R0 = state.render, M0 = state.mpr;
const g3 = $('#g-3d'), gM = $('#g-mpr'), gA = $('#g-arch'), gR = $('#g-recon');
const rebuild = () => requestBuild(false);
const refit = () => { if (state.src) state.src.fitTouched = true; rebuild(); };
const on3D = () => { syncPresetSeg(); moving(); request('gl'); };

addSelect(g3, 'r-mode', 'Render mode', R0, 'mode', [['surface', 'Surface (realistic)'], ['xray', 'X-ray glass (see-through)']], on3D);
addCheck(g3, 'r-gums', 'Show gums', R0, 'gums', on3D);
addRange(g3, 'r-guma', 'Gum opacity', R0, 'gumA', 0.1, 1, 0.01, v => v.toFixed(2), on3D);
addRange(g3, 'r-bonea', 'Bone opacity', R0, 'boneA', 0, 1, 0.01, v => v.toFixed(2), on3D);
addRange(g3, 'r-thr', 'Surface threshold', R0, 'thr', 0.15, 0.7, 0.005, v => v.toFixed(2), () => { requestMesh(); capKey = ''; on3D(); });
addSelect(g3, 'r-clip', 'Cut away at plane', R0, 'clip', [['none', 'None'], ['axial', 'Axial'], ['coronal', 'Coronal'], ['sagittal', 'Sagittal'], ['oblique', 'Cross-section / oblique']], on3D);
addCheck(g3, 'r-planes', 'Show slice planes', R0, 'planes', on3D);
addRow(g3, '<p class="note">Lower the bone opacity to see roots and implants. With a cut plane, the half facing you is removed.</p>');

const WPRESETS = { full: [255, 128], bone: [190, 150], teeth: [120, 190], soft: [110, 60] };
addSelect(gM, 's-wp', 'Window preset', M0, 'wpreset', [['full', 'Full range'], ['bone', 'Bone'], ['teeth', 'Teeth and enamel'], ['soft', 'Soft tissue'], ['custom', 'Custom']], () => { const w = WPRESETS[M0.wpreset]; if (w) { M0.win = w[0]; M0.lev = w[1]; syncControls(); } request('panes'); });
addRange(gM, 's-win', 'Window', M0, 'win', 10, 255, 1, v => String(v), () => { M0.wpreset = 'custom'; syncControls(); request('panes'); });
addRange(gM, 's-lev', 'Level', M0, 'lev', 0, 255, 1, v => String(v), () => { M0.wpreset = 'custom'; syncControls(); request('panes'); });
addSelect(gM, 's-obl', 'Fourth pane', M0, 'obl', [['arch', 'Cross-section across the arch'], ['free', 'Free oblique angle']], () => { updateOblUI(); request('gl', 'panes'); });
const yawW = addRange(gM, 's-yaw', 'Yaw', M0, 'yaw', -90, 90, 1, v => v + '°', () => request('gl', 'panes'));
const pitchW = addRange(gM, 's-pitch', 'Pitch', M0, 'pitch', -80, 80, 1, v => v + '°', () => request('gl', 'panes'));
function updateOblUI() { const free = M0.obl === 'free'; yawW.hidden = !free; pitchW.hidden = !free; $('#ob-title').textContent = free ? 'Oblique' : 'Cross-section'; }
updateOblUI();
{
  const row = addRow(gM, '<button class="btn small" type="button" id="btn-measure" aria-pressed="false">Measure distance</button><button class="btn small" type="button" id="btn-clear">Clear measurements</button>');
  row.querySelector('#btn-measure').addEventListener('click', () => setMeasure(!state.measure));
  row.querySelector('#btn-clear').addEventListener('click', () => { state.measures = []; request('panes'); });
  addRow(gM, '<p class="note">With Measure on, drag across a slice. Distances are approximate because OPG magnification varies.</p>');
}

addRange(gA, 'p-mid', 'Midline', P0, 'mid', 0.35, 0.65, 0.001, pct, refit);
addRange(gA, 'p-occl', 'Occlusal plane at midline', P0, 'occl', 0.3, 0.7, 0.001, v => pct(v) + ' down', refit);
addRange(gA, 'p-curve', 'Occlusal curve (smile +, frown \u2212)', P0, 'curve', -0.2, 0.4, 0.001, v => (v < 0 ? '\u2212' : '') + Math.abs(v).toFixed(3), refit);
addRange(gA, 'p-span', 'Condyle-to-condyle span', P0, 'span', 0.6, 1.3, 0.005, v => (v * 100).toFixed(0) + ' % of width', refit);
addRange(gA, 'p-archw', 'Arch width', P0, 'archW', 0.8, 1.25, 0.01, times, rebuild);
addRange(gA, 'p-archd', 'Arch depth', P0, 'archD', 0.8, 1.25, 0.01, times, rebuild);
addRange(gA, 'p-vscale', 'Vertical scale', P0, 'vscale', 0.7, 1.4, 0.01, times, rebuild);
{
  const row = addRow(gA, '<p class="note">Auto-fit sets the midline and occlusal curve. Set the span so the dashed blue lines sit on the condyles.</p><button class="btn small" type="button">Reset fit</button>');
  row.querySelector('button').addEventListener('click', () => { for (const k of ['mid', 'span', 'archW', 'archD', 'vscale']) P0[k] = DEFAULT_PARAMS[k]; syncControls(); rebuild(); });
}
addRange(gR, 'p-thick', 'Jaw thickness', P0, 'thick', 0.6, 1.5, 0.01, times, rebuild);
addRange(gR, 'p-round', 'Tooth roundness', P0, 'round', 0.5, 1.4, 0.01, v => v.toFixed(2), rebuild);
addRange(gR, 'p-floor', 'Soft-tissue cut-off', P0, 'floor', 0.05, 0.7, 0.005, v => v.toFixed(2), rebuild);
addRange(gR, 'p-shell', 'Enamel and cortex shell', P0, 'shell', 0, 1, 0.01, v => v.toFixed(2), rebuild);
addRange(gR, 'p-smooth', 'Smoothing', P0, 'smooth', 0, 4, 0.1, v => v.toFixed(1) + ' px', rebuild);
addSelect(gR, 'p-res', 'Resolution', P0, 'res', [['lite', 'Light (128 voxels across, phones)'], ['std', 'Standard (176 voxels across)'], ['high', 'High (240 voxels across, desktop)']], rebuild);
addCheck(gR, 'p-invert', 'Invert image (teeth look dark)', P0, 'invert', () => {
  if (!state.src) return;
  const S = state.src;
  prepSource(); S.sid = ++sidSeq; S.teeth = null; S.after = null;
  rebuild();
  if (wantModel()) { S.after = () => { if (state.src === S && S.teeth) { request('opg'); rebuild(); } }; requestSegment(); }
});
if (NETMETA) {
  const d = NETMETA.val_dice;
  addSelect(gR, 'p-teeth', 'Teeth detection', P0, 'teethSrc', [['model', 'Trained U-Net (real OPGs)'], ['contrast', 'Brightness contrast (no model)']], () => {
    const S = state.src; if (!S) return;
    if (P0.teethSrc === 'model' && !S.teeth) { if (!S.teethPending) { S.after = () => { request('opg'); rebuild(); }; requestSegment(); } }
    else { request('opg'); rebuild(); }
  });
  addCheck(gR, 'u-teeth', 'Tint detected teeth on the OPG', state.ui, 'showTeeth', () => request('opg'));
  addRow(gR, `<p class="note">The U-Net was trained on ${NETMETA.n_train || 96} real OPGs with tooth masks${d ? ` and scores a Dice of ${d.toFixed(2)} on ${NETMETA.n_val || 20} OPGs it never saw` : ''}. It finds teeth only; everything else is still estimated.</p>`);
}

/* =====================================================================
   Boot
   ===================================================================== */
buildChart();
// 2D canvases can be wiped by a GPU reset; redraw them when the browser restores them
for (const cv of [opgCv, glov, ...panes.map(p => p.cv)]) cv.addEventListener('contextrestored', () => { panes.forEach(p => { p.img = null; }); request('panes', 'opg', 'gl'); });
document.addEventListener('visibilitychange', () => { if (!document.hidden) request('gl', 'panes', 'opg'); });
initGL();
initWorker();
const ro = new ResizeObserver(scheduleResize);
ro.observe(views); document.querySelectorAll('.view').forEach(v => ro.observe(v)); ro.observe($('#drop'));
resizeAll();
loadDemo();
})();
