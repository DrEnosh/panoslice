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
    const target = P.res === 'high' ? 240 : 176;
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
  return { clamp, lerp, smoothstep, UPPER, LOWER, REF_L, thickAt, buildArch, nearestIdx, idxAtS, tableAt, occFracAt, gaussBlur, blurXY, percentileNormalize, buildVolume, makeNet: typeof makeNet === 'function' ? makeNet : null, segmentTeeth: typeof segmentTeeth === 'function' ? segmentTeeth : null };
}

function workerMain() {
  const SH = shared(); let src = null, net = null; const cache = {};
  const segment = sid => {
    if (!net || !src) { self.postMessage({ type: 'teeth', sid, error: 'no model' }); return; }
    try {
      const t0 = performance.now(), t = SH.segmentTeeth(net, src.norm, src.w, src.h);
      src.teeth = t;
      self.postMessage({ type: 'teeth', sid, prob: t.prob, w: t.w, h: t.h, ms: performance.now() - t0 });
    } catch (err) { self.postMessage({ type: 'teeth', sid, error: String((err && err.message) || err) }); }
  };
  self.onmessage = e => {
    const m = e.data;
    if (m.type === 'net') { net = SH.makeNet(m.meta, m.u16); return; }
    if (m.type === 'source') { src = m; src.teeth = null; cache.key = null; if (m.segment) segment(m.sid); return; }
    if (m.type === 'segment') { segment(m.sid); return; }
    if (m.type === 'build') {
      try {
        const r = SH.buildVolume(src.norm, src.raw, src.w, src.h, m.params, cache, m.params.teethSrc === 'model' ? src.teeth : null);
        self.postMessage({ type: 'built', id: m.id, r }, [r.rgba.buffer, r.slice.buffer]);
      } catch (err) { self.postMessage({ type: 'error', id: m.id, msg: String((err && err.message) || err) }); }
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
const DEFAULT_PARAMS = { mid: 0.5, span: 1, archW: 1, archD: 1, vscale: 1, occl: 0.515, curve: 0.075, thick: 1, floor: 0.28, shell: 0.5, smooth: 1, round: 0.85, res: 'std', invert: false, teethSrc: NETMETA ? 'model' : 'contrast' };
const PRESETS = {
  anatomy: { mode: 'surface', gums: true, boneA: 1 },
  bone: { mode: 'surface', gums: false, boneA: 1 },
  roots: { mode: 'surface', gums: false, boneA: 0.16 },
  xray: { mode: 'xray' },
};
const state = {
  src: null,
  params: { ...DEFAULT_PARAMS },
  render: { preset: 'anatomy', mode: 'surface', gums: true, gumA: 0.96, boneA: 1, thr: 0.33, opacity: 0.6, hq: true, clip: 'none', planes: false, spin: false },
  mpr: { win: 255, lev: 128, wpreset: 'full', obl: 'arch', yaw: 35, pitch: 20 },
  arch: null, vol: null,
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
let worker = null, buildSeq = 0, pending = null, current = null;
const mainCache = {};
function initWorker() {
  try {
    const src = `const shared = ${shared.toString()};\n(${workerMain.toString()})();`;
    worker = new Worker(URL.createObjectURL(new Blob([src], { type: 'text/javascript' })));
    netReady.then(ok => { if (ok && worker) worker.postMessage({ type: 'net', meta: NETMETA, u16: NETU16.slice() }); });
    worker.onmessage = e => {
      const m = e.data;
      if (m.type === 'teeth') { onTeeth(m); return; }
      if (!current || m.id !== current.id) return;
      if (m.type === 'built') finishBuild(m.r, current);
      else { console.error(m.msg); worker = null; runJob(current); }
    };
    worker.onerror = e => {
      console.warn('Worker unavailable, working on the main thread.', e.message || ''); e.preventDefault && e.preventDefault(); worker = null;
      if (state.src && state.src.teethPending) netReady.then(segmentOnMain); else if (current) runJob(current);
    };
  } catch (e) { worker = null; }
}
function wantModel() { return !!NETMETA && state.params.teethSrc === 'model'; }
function sendSource(segment) {
  const S = state.src; if (!S) return;
  const post = seg => {
    if (worker) worker.postMessage({ type: 'source', sid: S.sid, segment: seg, norm: S.norm, raw: S.prep, w: S.w, h: S.h });
    else if (seg) setTimeout(segmentOnMain, 30);
  };
  if (!segment) { post(false); return; }
  S.teethPending = true; S.teeth = null;
  netReady.then(ok => {                       // the weights may still be downloading
    if (state.src !== S) return;
    if (ok) post(true); else { post(false); onTeeth({ sid: S.sid, error: 'no model' }); }
  });
}
function requestSegment() {
  const S = state.src; if (!S) return;
  S.teethPending = true;
  netReady.then(ok => {
    if (state.src !== S) return;
    if (!ok) { onTeeth({ sid: S.sid, error: 'no model' }); return; }
    if (worker) worker.postMessage({ type: 'segment', sid: S.sid }); else setTimeout(segmentOnMain, 30);
  });
}
function segmentOnMain() {
  const S = state.src, net = getMainNet(); if (!S) return;
  if (!net) { onTeeth({ sid: S.sid, error: 'no model' }); return; }
  const t0 = performance.now();
  try { const t = SH.segmentTeeth(net, S.norm, S.w, S.h); onTeeth({ sid: S.sid, prob: t.prob, w: t.w, h: t.h, ms: performance.now() - t0 }); }
  catch (e) { console.error(e); onTeeth({ sid: S.sid, error: String(e) }); }
}
/* Teeth probability arrived: tint the thumbnail, fit the arch from it, then build. */
function onTeeth(m) {
  const S = state.src; if (!S || m.sid !== S.sid) return;
  S.teethPending = false;
  if (m.error) { S.teeth = null; state.params.teethSrc = 'contrast'; syncControls(); state.notice = 'The teeth model couldn\u2019t run here, so teeth were separated by brightness instead.'; }
  else { S.teeth = { prob: m.prob, w: m.w, h: m.h }; S.teethMs = m.ms; buildTeethOverlay(); }
  if (S.after) { const f = S.after; S.after = null; f(); }
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
  pending = { params: { ...state.params }, reset: !!reset || !!(pending && pending.reset) };
  if (!current) nextJob();
}
function nextJob() {
  const job = pending; pending = null;
  if (!job || !state.src) { showBusy(false); return; }
  job.id = ++buildSeq; current = job; showBusy(true); setStatus('Building the 3D model…');
  runJob(job);
}
function runJob(job) {
  if (worker) { worker.postMessage({ type: 'build', id: job.id, params: job.params }); return; }
  setTimeout(() => {
    if (current !== job) return;
    try { finishBuild(SH.buildVolume(state.src.norm, state.src.prep, state.src.w, state.src.h, job.params, mainCache, job.params.teethSrc === 'model' ? state.src.teeth : null), job); }
    catch (e) { console.error(e); setStatus('The reconstruction failed on this image. Try another export of the OPG.', true); current = null; showBusy(false); }
  }, 16);
}
function finishBuild(r, job) {
  current = null;
  state.vol = r;
  state.arch = SH.buildArch(r.archW, r.archD);
  uploadVolume();
  if (job.reset || !state.P) defaultCrosshair(); else clampP();
  const tsrc = job.params.teethSrc === 'model' && state.src.teeth ? ` · teeth by U-Net (${(state.src.teethMs / 1000).toFixed(1)} s)` : '';
  setStatus(`${r.nx}×${r.ny}×${r.nz} voxels · ${r.vox.toFixed(2)} mm · ${Math.round(r.ms)} ms${tsrc}`);
  if (state.notice) { setStatus(state.notice, true); state.notice = null; }
  request('gl', 'panes', 'opg', 'readout');
  if (pending) nextJob(); else showBusy(false);
}
function showBusy(on) { $('#busy').hidden = !on; }
function setStatus(msg, warn) { const el = $('#status'); el.textContent = msg; el.title = msg; el.classList.toggle('warn', !!warn); }

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
    sendSource(false); requestBuild(true);
    if (wantModel()) { state.src.after = () => { if (state.params.teethSrc === 'model') requestBuild(false); }; requestSegment(); }
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
  const sc = Math.min(1, 1600 / img.naturalWidth, 1000 / img.naturalHeight);
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
  const aspect = w / h;
  const finish = () => {
    const fit = fitArch(true);
    syncControls();
    if (aspect < 1.45 || w < 500) state.notice = `This image is ${w}\u00d7${h}. OPGs are usually about twice as wide as tall, so the 3D model may not make sense.`;
    else if (fit && !fit.ok) state.notice = 'Auto-fit wasn\u2019t confident on this image. Check that the tooth ticks sit on the teeth and adjust Arch fit if they don\u2019t.';
    requestBuild(true);
  };
  if (wantModel()) { setStatus('Finding teeth with the trained model\u2026'); showBusy(true); state.src.after = finish; sendSource(true); }
  else { sendSource(false); finish(); }
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
   WebGL2 renderer
   ===================================================================== */
const glc = $('#gl'), glov = $('#gl-ov'), govx = glov.getContext('2d');
const gl = glc.getContext('webgl2', { alpha: false, antialias: false, preserveDrawingBuffer: false });
let prog, lineProg, vaoQuad, vaoLine, lineBuf, tex, U = {}, UL = {}, glOK = !!gl;
const VS_QUAD = `#version 300 es
layout(location=0) in vec2 aPos; out vec2 vNdc;
void main(){ vNdc = aPos; gl_Position = vec4(aPos, 0.0, 1.0); }`;
const FS_VOL = `#version 300 es
precision highp float;
precision highp sampler3D;
in vec2 vNdc; out vec4 o;
uniform sampler3D uVol;
uniform mat4 uInvVP;
uniform vec3 uBMin, uBMax, uTexel, uKey, uClipN;
uniform int uMode, uClipOn, uGums, uHQ;
uniform float uThr, uStep, uBoneA, uGumA, uOpacity, uClipD;

vec4 T(vec3 p){ return texture(uVol, (p - uBMin) / (uBMax - uBMin)); }
vec3 grad(vec3 p, int c){
  vec3 h = (uBMax - uBMin) * uTexel * 1.25;
  return vec3(T(p + vec3(h.x,0.,0.))[c] - T(p - vec3(h.x,0.,0.))[c],
              T(p + vec3(0.,h.y,0.))[c] - T(p - vec3(0.,h.y,0.))[c],
              T(p + vec3(0.,0.,h.z))[c] - T(p - vec3(0.,0.,h.z))[c]) / (2.0 * h);
}
vec3 normalAt(vec3 p, int c, vec3 rd){ vec3 n = -normalize(grad(p, c) + vec3(1e-6)); return dot(n, rd) > 0.0 ? -n : n; }
float refine(vec3 ro, vec3 rd, float a, float b, int c, float th){
  for (int k = 0; k < 6; k++){ float m = 0.5 * (a + b); if (T(ro + rd * m)[c] >= th) b = m; else a = m; }
  return b;
}
float hash13(vec3 p){ p = fract(p * 0.1031); p += dot(p, p.zyx + 31.32); return fract((p.x + p.y) * p.z); }
float vnoise(vec3 p){
  vec3 i = floor(p), f = fract(p); f = f * f * (3.0 - 2.0 * f);
  return mix(mix(mix(hash13(i), hash13(i + vec3(1,0,0)), f.x), mix(hash13(i + vec3(0,1,0)), hash13(i + vec3(1,1,0)), f.x), f.y),
             mix(mix(hash13(i + vec3(0,0,1)), hash13(i + vec3(1,0,1)), f.x), mix(hash13(i + vec3(0,1,1)), hash13(i + vec3(1,1,1)), f.x), f.y), f.z);
}
vec3 sky(vec3 d){ float t = clamp(d.z * 0.5 + 0.5, 0.0, 1.0); return mix(vec3(0.05,0.045,0.04), vec3(0.85,0.9,1.0), pow(t, 1.4)); }
float occ(vec3 p){
  vec4 s = T(p); float hard = smoothstep(uThr - 0.05, uThr + 0.15, s.r);
  return max(hard * smoothstep(0.35, 0.65, max(s.g, s.a)), hard * uBoneA);
}
const vec3 AOK[12] = vec3[12](
  vec3(0.,0.526,0.851), vec3(0.,-0.526,0.851), vec3(0.,0.526,-0.851), vec3(0.,-0.526,-0.851),
  vec3(0.526,0.851,0.), vec3(-0.526,0.851,0.), vec3(0.526,-0.851,0.), vec3(-0.526,-0.851,0.),
  vec3(0.851,0.,0.526), vec3(-0.851,0.,0.526), vec3(0.851,0.,-0.526), vec3(-0.851,0.,-0.526));
float ambOcc(vec3 p, vec3 n){
  float o = 0.0, tot = 0.0; int ns = uHQ == 1 ? 12 : 6;
  for (int i = 0; i < 12; i++){
    if (i >= ns) break;
    vec3 d = normalize(n + 0.85 * AOK[i]);
    o += occ(p + d * 1.1) * 0.55 + occ(p + d * 2.8) * 0.45; tot += 1.0;
  }
  return clamp(1.0 - 1.15 * o / tot, 0.12, 1.0);
}
float shadow(vec3 p){
  if (uHQ == 0) return 1.0;
  float res = 1.0, t = 0.9;
  for (int i = 0; i < 40; i++){ if (t > 34.0 || res < 0.04) break; res *= 1.0 - 0.5 * occ(p + uKey * t); t += 0.55 + t * 0.05; }
  return res;
}
vec3 light(vec3 p, vec3 n, vec3 V, vec3 alb, float specI, float specP, float wrap, float metal){
  float ao = ambOcc(p, n), sh = shadow(p + n * 0.6);
  float ndl = dot(n, uKey);
  float diff = max((ndl + wrap) / (1.0 + wrap), 0.0);
  vec3 L2 = normalize(vec3(-uKey.x, -uKey.y, 0.4));
  float fill = max(dot(n, L2), 0.0);
  float nv = max(dot(n, V), 0.0);
  vec3 hemi = mix(vec3(0.10,0.09,0.085), vec3(0.34,0.37,0.42), n.z * 0.5 + 0.5);
  vec3 H = normalize(uKey + V);
  float spec = pow(max(dot(n, H), 0.0), specP) * specI * sh * (0.4 + 0.6 * smoothstep(0.0, 0.2, ndl));
  float fres = 0.04 + 0.96 * pow(1.0 - nv, 5.0);
  vec3 env = sky(reflect(-V, n));
  vec3 c = alb * (hemi * ao + 2.3 * diff * sh * vec3(1.0,0.96,0.9) + 0.4 * fill * ao * vec3(0.72,0.8,1.0));
  vec3 m = alb * (env * 1.3 * ao + 0.6 * diff * sh) + spec * 2.0;
  c = mix(c + spec * vec3(1.0,0.98,0.94) + fres * env * 0.22 * ao, m, metal);
  return c + pow(1.0 - nv, 3.0) * 0.12 * ao * vec3(0.55,0.65,0.8);
}
vec3 aces(vec3 x){ return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0); }
vec3 enamelC(float v){ return mix(vec3(0.8,0.7,0.52), vec3(0.9,0.86,0.76), smoothstep(0.5, 0.8, v)); }
const vec3 BONE = vec3(0.72,0.6,0.43);
const vec3 GUM = vec3(0.46,0.1,0.11);
const vec3 METAL = vec3(0.62,0.64,0.68);
vec3 shadeHard(vec3 p, vec3 n, vec3 V, int mat, float dens, bool root){
  if (mat == 2) return light(p, n, V, METAL, 1.0, 90.0, 0.0, 1.0);
  if (mat == 1) return light(p, n, V, root ? vec3(0.82,0.7,0.52) : enamelC(dens), root ? 0.18 : 0.55, root ? 30.0 : 70.0, 0.35, 0.0);
  float nz = vnoise(p * 1.6) * 0.6 + vnoise(p * 4.1) * 0.4;
  return light(p, n, V, BONE * (0.8 + 0.28 * nz), 0.06, 16.0, 0.1, 0.0);
}
vec3 bgCol(){
  vec2 q = vNdc; float r = length(q * vec2(0.9, 1.0));
  vec3 top = vec3(0.105,0.13,0.15), bot = vec3(0.035,0.045,0.055);
  return mix(mix(bot, top, clamp(q.y * 0.5 + 0.6, 0.0, 1.0)), bot * 0.8, smoothstep(0.6, 1.5, r));
}
vec3 toSRGB(vec3 c){ return pow(aces(c), vec3(1.0 / 2.2)); }

void main(){
  vec3 BG = bgCol();
  vec4 a = uInvVP * vec4(vNdc, -1.0, 1.0); a /= a.w;
  vec4 b = uInvVP * vec4(vNdc, 1.0, 1.0); b /= b.w;
  vec3 ro = a.xyz, rd = normalize(b.xyz - a.xyz), V = -rd;
  vec3 inv = 1.0 / rd;
  vec3 t0 = (uBMin - ro) * inv, t1 = (uBMax - ro) * inv;
  vec3 lo3 = min(t0, t1), hi3 = max(t0, t1);
  float tn = max(max(max(lo3.x, lo3.y), lo3.z), 0.0);
  float tf = min(min(hi3.x, hi3.y), hi3.z);
  vec3 bgOut = uMode == 2 ? vec3(0.02,0.025,0.03) : BG;
  if (tf <= tn) { o = vec4(bgOut, 1.0); return; }
  bool cap = false;
  if (uClipOn == 1) {
    float dn = dot(uClipN, rd), dr = dot(uClipN, ro) - uClipD;
    if (abs(dn) < 1e-6) { if (dr > 0.0) { o = vec4(bgOut, 1.0); return; } }
    else {
      float tc = -dr / dn;
      if (dn > 0.0) tf = min(tf, tc);
      else if (tc > tn) { tn = tc; cap = true; }
    }
    if (tf <= tn) { o = vec4(bgOut, 1.0); return; }
  }
  if (uMode == 0) {
    vec3 acc = vec3(0.0); float A = 0.0;
    bool inHard = false, inGum = false, first = true;
    float t = tn, prev = tn, prevG = 0.0;
    for (int i = 0; i < 2400; i++) {
      if (t > tf || A > 0.985) break;
      vec3 p = ro + rd * t; vec4 s = T(p);
      bool hard = s.r >= uThr;
      if (first && cap) {
        if (hard) {
          int mat = s.a > 0.5 ? 2 : (s.g >= 0.5 ? 1 : 0);
          vec3 c = (mat == 2 ? METAL * 1.3 : (mat == 1 ? enamelC(s.r) : BONE)) * (0.3 + 0.8 * s.r) * 1.3;
          acc += (1.0 - A) * c; A = 1.0; break;
        }
        if (uGums == 1 && s.b >= 0.5) { float ga = uGumA * 0.85; acc += (1.0 - A) * ga * GUM * 1.2; A += (1.0 - A) * ga; inGum = true; }
      } else if (uGums == 1) {
        bool g = s.b >= 0.5;
        if (g && !inGum) {
          float th = refine(ro, rd, prev, t, 2, 0.5); vec3 q = ro + rd * th;
          vec3 n = normalAt(q, 2, rd);
          float nv = max(dot(n, V), 0.0);
          vec3 c = light(q, n, V, GUM, 0.22, 70.0, 0.5, 0.0) + vec3(0.5,0.08,0.08) * 0.35 * pow(1.0 - nv, 2.0);
          float ga = clamp(uGumA * (0.85 + 0.15 * (1.0 - nv)), 0.0, 1.0);
          acc += (1.0 - A) * ga * c; A += (1.0 - A) * ga;
        }
        inGum = g;
      }
      if (hard && !inHard) {
        float th = refine(ro, rd, prev, t, 0, uThr); vec3 q = ro + rd * th;
        vec4 sq = T(q + rd * 0.45);
        int mat = sq.a > 0.5 ? 2 : (sq.g >= 0.5 ? 1 : 0);
        vec3 n = normalAt(q, 0, rd);
        if (mat != 0) { acc += (1.0 - A) * shadeHard(q, n, V, mat, sq.r, false); A = 1.0; break; }
        if (uBoneA > 0.003) { acc += (1.0 - A) * uBoneA * shadeHard(q, n, V, 0, sq.r, false); A += (1.0 - A) * uBoneA; }
        inHard = true;
      } else if (hard && inHard && max(s.g, s.a) >= 0.5 && prevG < 0.5) {
        float th = refine(ro, rd, prev, t, 1, 0.5); vec3 q = ro + rd * th;
        vec4 sq = T(q + rd * 0.3);
        acc += (1.0 - A) * shadeHard(q, normalAt(q, 1, rd), V, sq.a > 0.5 ? 2 : 1, sq.r, true); A = 1.0; break;
      } else if (!hard) inHard = false;
      first = false; prev = t; prevG = max(s.g, s.a); t += uStep;
    }
    vec3 col = A > 0.0 ? toSRGB(acc / A) : vec3(0.0);
    o = vec4(mix(BG, col, A), 1.0); return;
  }
  float t = tn + hash13(vec3(gl_FragCoord.xy, 1.0)) * uStep;
  if (uMode == 1) {
    vec3 acc = vec3(0.0); float A = 0.0;
    for (int i = 0; i < 2400; i++) {
      if (t > tf || A > 0.98) break;
      vec3 p = ro + rd * t; vec4 s = T(p);
      float hard = smoothstep(uThr - 0.04, uThr + 0.2, s.r);
      float th = smoothstep(0.4, 0.6, max(s.g, s.a)) * hard;
      float bo = hard * (1.0 - th) * uBoneA;
      float gu = uGums == 1 ? smoothstep(0.4, 0.6, s.b) * (1.0 - hard) * uGumA : 0.0;
      float al = (th * 0.85 + bo * 0.45 + gu * 0.12) * uOpacity;
      if (al > 0.002) {
        al = 1.0 - pow(1.0 - clamp(al, 0.0, 0.995), uStep / 0.5);
        vec3 alb = (th * (s.a > 0.5 ? METAL : enamelC(s.r)) + bo * BONE + gu * GUM) / max(th + bo + gu, 1e-4);
        vec3 c = alb * (0.3 + 0.9 * s.r);
        if (uHQ == 1 && hard > 0.05) {
          vec3 g = grad(p, 0); float gm = length(g);
          if (gm > 0.02) { vec3 n = -g / gm; if (dot(n, rd) > 0.0) n = -n; c = alb * (0.3 + 1.6 * max(dot(n, uKey), 0.0)) + pow(max(dot(n, normalize(uKey + V)), 0.0), 40.0) * 0.25; }
        }
        acc += (1.0 - A) * al * c; A += (1.0 - A) * al;
      }
      t += uStep;
    }
    vec3 col = A > 0.0 ? toSRGB(acc / A) : vec3(0.0);
    o = vec4(mix(BG, col, A), 1.0); return;
  }
  float sum = 0.0;
  for (int i = 0; i < 2400; i++) { if (t > tf) break; vec4 s = T(ro + rd * t); sum += (s.r + 0.12 * s.b) * uStep; t += uStep; }
  float g = smoothstep(0.12, 0.8, 1.0 - exp(-sum * 0.11));
  o = vec4(mix(vec3(0.02,0.025,0.03), vec3(0.93,0.95,0.97), g), 1.0);
}`;
const VS_LINE = `#version 300 es
layout(location=0) in vec3 aPos; uniform mat4 uVP;
void main(){ gl_Position = uVP * vec4(aPos, 1.0); }`;
const FS_LINE = `#version 300 es
precision mediump float; uniform vec4 uCol; out vec4 o;
void main(){ o = uCol; }`;

function compile(vs, fs) {
  const mk = (type, src) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s)); return s; };
  const p = gl.createProgram(); gl.attachShader(p, mk(gl.VERTEX_SHADER, vs)); gl.attachShader(p, mk(gl.FRAGMENT_SHADER, fs)); gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
function initGL() {
  if (!gl) { showGLMessage('The 3D view needs WebGL 2, which this browser doesn’t offer. The slice views still work.'); return; }
  try { prog = compile(VS_QUAD, FS_VOL); lineProg = compile(VS_LINE, FS_LINE); }
  catch (e) { console.error(e); glOK = false; showGLMessage('The 3D renderer couldn’t start on this device. The slice views still work.'); return; }
  for (const n of ['uVol', 'uInvVP', 'uBMin', 'uBMax', 'uTexel', 'uKey', 'uClipN', 'uMode', 'uClipOn', 'uGums', 'uHQ', 'uThr', 'uStep', 'uBoneA', 'uGumA', 'uOpacity', 'uClipD']) U[n] = gl.getUniformLocation(prog, n);
  for (const n of ['uVP', 'uCol']) UL[n] = gl.getUniformLocation(lineProg, n);
  vaoQuad = gl.createVertexArray(); gl.bindVertexArray(vaoQuad);
  const qb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, qb);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
  vaoLine = gl.createVertexArray(); gl.bindVertexArray(vaoLine);
  lineBuf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf);
  gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
  gl.bindVertexArray(null);
  tex = gl.createTexture();
  glc.addEventListener('webglcontextlost', e => { e.preventDefault(); glOK = false; showGLMessage('The 3D view lost its graphics context. Reload the page to restore it.'); });
}
function showGLMessage(msg) { const el = $('#gl-msg'); el.textContent = msg; el.hidden = false; }
function uploadVolume() {
  if (!glOK) return;
  const V = state.vol;
  gl.bindTexture(gl.TEXTURE_3D, tex);
  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  gl.texImage3D(gl.TEXTURE_3D, 0, gl.RGBA8, V.nx, V.ny, V.nz, 0, gl.RGBA, gl.UNSIGNED_BYTE, V.rgba);
  gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  for (const w of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T, gl.TEXTURE_WRAP_R]) gl.texParameteri(gl.TEXTURE_3D, w, gl.CLAMP_TO_EDGE);
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
  const vp = mul4(perspective(FOV, aspect, dist * 0.05, dist * 3), view);
  const right = [view[0], view[4], view[8]], up = [view[1], view[5], view[9]], back = [view[2], view[6], view[10]];
  const key = V3.norm(V3.add(V3.add(V3.mul(right, -0.5), V3.mul(up, 0.75)), V3.mul(back, 0.55)));
  return { eye, view, vp, inv: invert4(vp), dist, diag, key, C };
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
  const sc = state.moving ? Math.min(dpr, 1) * 0.7 : Math.min(dpr, 2, Math.sqrt(2.6e6 / area));
  const w = Math.max(1, Math.round(glc.clientWidth * sc)), h = Math.max(1, Math.round(glc.clientHeight * sc));
  if (glc.width !== w || glc.height !== h) { glc.width = w; glc.height = h; }
}
function render3D() {
  if (!glOK || !state.vol || !glc.clientWidth) return;
  sizeGL();
  const w = glc.width, h = glc.height;
  const V = state.vol, M = camera(w / h), R = state.render;
  gl.viewport(0, 0, w, h);
  gl.disable(gl.BLEND);
  gl.useProgram(prog);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_3D, tex);
  gl.uniform1i(U.uVol, 0);
  gl.uniformMatrix4fv(U.uInvVP, false, M.inv);
  gl.uniform3fv(U.uBMin, V.bmin); gl.uniform3fv(U.uBMax, V.bmax);
  gl.uniform3f(U.uTexel, 1 / V.nx, 1 / V.ny, 1 / V.nz);
  gl.uniform3fv(U.uKey, M.key);
  gl.uniform1i(U.uMode, R.mode === 'surface' ? 0 : R.mode === 'volume' ? 1 : 2);
  gl.uniform1i(U.uGums, R.gums ? 1 : 0);
  const hq = R.hq && !state.moving;
  gl.uniform1i(U.uHQ, hq ? 1 : 0);
  gl.uniform1f(U.uThr, R.thr); gl.uniform1f(U.uOpacity, R.opacity);
  gl.uniform1f(U.uBoneA, R.boneA); gl.uniform1f(U.uGumA, R.gumA);
  gl.uniform1f(U.uStep, V.vox * (state.moving ? 1.1 : 0.45));
  const cp = clipPlane(M);
  gl.uniform1i(U.uClipOn, cp ? 1 : 0);
  gl.uniform3fv(U.uClipN, cp ? cp.n : [0, 0, 1]); gl.uniform1f(U.uClipD, cp ? cp.d : 0);
  gl.bindVertexArray(vaoQuad); gl.drawArrays(gl.TRIANGLES, 0, 3);
  drawLines3D(M);
  gl.bindVertexArray(null);
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
  gl.useProgram(lineProg); gl.bindVertexArray(vaoLine);
  gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(verts), gl.DYNAMIC_DRAW);
  gl.uniformMatrix4fv(UL.uVP, false, M.vp);
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  for (const g of groups) {
    const c = PLANE_RGB[g.key];
    if (g.kind === 'fill') { gl.uniform4f(UL.uCol, c[0], c[1], c[2], 0.035); gl.drawArrays(gl.TRIANGLES, g.start, g.count); }
    else { gl.uniform4f(UL.uCol, c[0], c[1], c[2], g.kind === 'mark' ? 1 : 0.8); gl.drawArrays(gl.LINES, g.start, g.count); }
  }
  gl.disable(gl.BLEND);
}
function drawTriad(M) {
  const w = glov.width, h = glov.height, d = w / Math.max(1, glov.clientWidth);
  govx.clearRect(0, 0, w, h);
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
    const dp = Math.min(window.devicePixelRatio || 1, Math.sqrt(260000 / Math.max(1, cw * chh)), 2);
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
const on3D = () => { syncPresetSeg(); moving(); request('gl'); };

addSelect(g3, 'r-mode', 'Render mode', R0, 'mode', [['surface', 'Surface (realistic)'], ['volume', 'Volume (see-through density)'], ['xray', 'X-ray (simulated radiograph)']], on3D);
addCheck(g3, 'r-gums', 'Show gums', R0, 'gums', on3D);
addRange(g3, 'r-guma', 'Gum opacity', R0, 'gumA', 0.1, 1, 0.01, v => v.toFixed(2), on3D);
addRange(g3, 'r-bonea', 'Bone opacity', R0, 'boneA', 0, 1, 0.01, v => v.toFixed(2), on3D);
addRange(g3, 'r-thr', 'Surface threshold', R0, 'thr', 0.1, 0.8, 0.005, v => v.toFixed(2), on3D);
addRange(g3, 'r-opac', 'Volume density', R0, 'opacity', 0.1, 1.5, 0.01, v => v.toFixed(2), on3D);
addCheck(g3, 'r-hq', 'Soft shadows and ambient occlusion', R0, 'hq', on3D);
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

addRange(gA, 'p-mid', 'Midline', P0, 'mid', 0.35, 0.65, 0.001, pct, rebuild);
addRange(gA, 'p-occl', 'Occlusal plane at midline', P0, 'occl', 0.3, 0.7, 0.001, v => pct(v) + ' down', rebuild);
addRange(gA, 'p-curve', 'Occlusal curve (smile +, frown \u2212)', P0, 'curve', -0.2, 0.4, 0.001, v => (v < 0 ? '\u2212' : '') + Math.abs(v).toFixed(3), rebuild);
addRange(gA, 'p-span', 'Condyle-to-condyle span', P0, 'span', 0.6, 1.3, 0.005, v => (v * 100).toFixed(0) + ' % of width', rebuild);
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
addSelect(gR, 'p-res', 'Resolution', P0, 'res', [['std', 'Standard (176 voxels across)'], ['high', 'High (240 voxels across, desktop)']], rebuild);
addCheck(gR, 'p-invert', 'Invert image (teeth look dark)', P0, 'invert', () => {
  if (!state.src) return;
  prepSource();
  if (wantModel()) { state.src.after = rebuild; setStatus('Finding teeth with the trained model\u2026'); showBusy(true); sendSource(true); } else { sendSource(false); rebuild(); }
});
if (NETMETA) {
  const d = NETMETA.val_dice;
  addSelect(gR, 'p-teeth', 'Teeth detection', P0, 'teethSrc', [['model', 'Trained U-Net (real OPGs)'], ['contrast', 'Brightness contrast (no model)']], () => {
    const S = state.src; if (!S) return;
    if (P0.teethSrc === 'model' && !S.teeth) { S.after = () => { request('opg'); rebuild(); }; setStatus('Finding teeth with the trained model\u2026'); showBusy(true); requestSegment(); }
    else { request('opg'); rebuild(); }
  });
  addCheck(gR, 'u-teeth', 'Tint detected teeth on the OPG', state.ui, 'showTeeth', () => request('opg'));
  addRow(gR, `<p class="note">The U-Net was trained on ${NETMETA.n_train || 96} real OPGs with tooth masks${d ? ` and scores a Dice of ${d.toFixed(2)} on ${NETMETA.n_val || 20} OPGs it never saw` : ''}. It finds teeth only; everything else is still estimated.</p>`);
}

/* =====================================================================
   Boot
   ===================================================================== */
buildChart();
initGL();
initWorker();
const ro = new ResizeObserver(scheduleResize);
ro.observe(views); document.querySelectorAll('.view').forEach(v => ro.observe(v)); ro.observe($('#drop'));
resizeAll();
loadDemo();
})();
