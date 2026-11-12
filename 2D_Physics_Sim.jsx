/* 2D Physics panel for EffectCraft: rigid-body physics baked to keyframes. See README.md */
/* ---------------------------------------------------------------------------
   PhysEngine: 2D rigid-body physics for convex polygons.
   SAT + clipping contact generation, sequential-impulse solver with warm
   starting, Coulomb friction, restitution. Pure JS (ES5), no dependencies.
   Coordinates: pixels, y down. Angles: radians, positive = clockwise on screen.
--------------------------------------------------------------------------- */
var PhysEngine = (function () {
  function cross2(ax, ay, bx, by) { return ax * by - ay * bx; }

  function hull(pts) {                       // convex hull, result has positive shoelace area
    var p = pts.slice().sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; });
    var q = [], i;
    for (i = 0; i < p.length; i++) {
      var l = q[q.length - 1];
      if (!l || Math.abs(p[i][0] - l[0]) > 1e-6 || Math.abs(p[i][1] - l[1]) > 1e-6) q.push(p[i]);
    }
    if (q.length < 3) return q;
    function cr(o, a, b) { return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); }
    var lo = [], up = [];
    for (i = 0; i < q.length; i++) {
      while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], q[i]) <= 0) lo.pop();
      lo.push(q[i]);
    }
    for (i = q.length - 1; i >= 0; i--) {
      while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], q[i]) <= 0) up.pop();
      up.push(q[i]);
    }
    lo.pop(); up.pop();
    return lo.concat(up);
  }

  function decimate(h, maxN) {               // keep convexity by sub-sampling hull vertices
    if (h.length <= maxN) return h;
    var out = [];
    for (var i = 0; i < maxN; i++) out.push(h[Math.floor(i * h.length / maxN)]);
    return out;
  }

  function polyProps(v) {                    // area, centroid, inertia about centroid (density 1)
    var a = 0, cx = 0, cy = 0, I = 0, n = v.length, i, j, cr;
    for (i = 0; i < n; i++) {
      j = (i + 1) % n;
      cr = cross2(v[i][0], v[i][1], v[j][0], v[j][1]);
      a += cr; cx += (v[i][0] + v[j][0]) * cr; cy += (v[i][1] + v[j][1]) * cr;
      I += cr * (v[i][0] * v[i][0] + v[i][0] * v[j][0] + v[j][0] * v[j][0] +
                 v[i][1] * v[i][1] + v[i][1] * v[j][1] + v[j][1] * v[j][1]) / 12;
    }
    a /= 2; cx /= (6 * a); cy /= (6 * a);
    var area = Math.abs(a);
    return { area: area, cx: cx, cy: cy, inertia: Math.max(1e-6, Math.abs(I) - area * (cx * cx + cy * cy)) };
  }


  function area2(p) { var a = 0, n = p.length, i, j; for (i = 0; i < n; i++) { j = (i + 1) % n; a += p[i][0] * p[j][1] - p[j][0] * p[i][1]; } return a / 2; }
  function cleanPoly(p) {                    // drop duplicate and collinear vertices
    var q = [], i, l;
    for (i = 0; i < p.length; i++) { l = q[q.length - 1]; if (!l || Math.abs(p[i][0] - l[0]) > 1e-4 || Math.abs(p[i][1] - l[1]) > 1e-4) q.push(p[i]); }
    while (q.length > 1 && Math.abs(q[0][0] - q[q.length - 1][0]) <= 1e-4 && Math.abs(q[0][1] - q[q.length - 1][1]) <= 1e-4) q.pop();
    var again = true;
    while (again && q.length > 3) {
      again = false;
      for (i = 0; i < q.length; i++) {
        var a = q[(i + q.length - 1) % q.length], b = q[i], c = q[(i + 1) % q.length];
        if (Math.abs((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) < 1e-3) { q.splice(i, 1); again = true; break; }
      }
    }
    return q;
  }
  function simplifyPoly(p, maxN) {           // remove the least significant vertices
    var q = p.slice();
    while (q.length > 3) {
      var best = -1, bestA = Infinity, n = q.length, i;
      for (i = 0; i < n; i++) {
        var a = q[(i + n - 1) % n], b = q[i], c = q[(i + 1) % n];
        var ar = Math.abs((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) / 2;
        if (ar < bestA) { bestA = ar; best = i; }
      }
      if (q.length <= maxN && bestA > 0.5) break;
      q.splice(best, 1);
    }
    return q;
  }
  function orient(a, b, p) { return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]); }
  function triangulate(p) {                  // ear clipping, p must be positive-area; returns index triples or null
    var idx = [], tris = [], i, guard = 0;
    for (i = 0; i < p.length; i++) idx.push(i);
    while (idx.length > 3 && guard++ < 20000) {
      var found = false, m = idx.length, k;
      for (k = 0; k < m && !found; k++) {
        var ia = idx[(k + m - 1) % m], ib = idx[k], ic = idx[(k + 1) % m], a = p[ia], b = p[ib], c = p[ic];
        if (orient(a, b, c) <= 1e-9) continue;
        var ear = true;
        for (var j = 0; j < m; j++) {
          var id = idx[j]; if (id === ia || id === ib || id === ic) continue;
          var q = p[id];
          if ((q[0] === a[0] && q[1] === a[1]) || (q[0] === b[0] && q[1] === b[1]) || (q[0] === c[0] && q[1] === c[1])) continue;
          if (orient(a, b, q) >= 0 && orient(b, c, q) >= 0 && orient(c, a, q) >= 0) { ear = false; break; }
        }
        if (ear) { tris.push([ia, ib, ic]); idx.splice(k, 1); found = true; }
      }
      if (!found) return null;
    }
    if (idx.length === 3) tris.push([idx[0], idx[1], idx[2]]);
    return tris;
  }
  function isConvexIdx(p, ids) {
    var n = ids.length, i;
    for (i = 0; i < n; i++) if (orient(p[ids[i]], p[ids[(i + 1) % n]], p[ids[(i + 2) % n]]) < -1e-7) return false;
    return true;
  }
  function mergePieces(p, tris) {            // Hertel-Mehlhorn: merge triangles while the result stays convex
    var pc = tris.slice(), changed = true, a, b, ia, ib;
    while (changed) {
      changed = false;
      for (a = 0; a < pc.length && !changed; a++) for (b = a + 1; b < pc.length && !changed; b++) {
        var A = pc[a], B = pc[b], la = A.length, lb = B.length;
        for (ia = 0; ia < la && !changed; ia++) {
          var u = A[ia], v = A[(ia + 1) % la];
          for (ib = 0; ib < lb; ib++) if (B[ib] === v && B[(ib + 1) % lb] === u) {
            var mg = [], k;
            for (k = 1; k <= la; k++) mg.push(A[(ia + k) % la]);
            for (k = 2; k < lb; k++) mg.push(B[(ib + k) % lb]);
            if (isConvexIdx(p, mg)) { pc[a] = mg; pc.splice(b, 1); changed = true; }
            break;
          }
        }
      }
    }
    return pc;
  }
  // concave polygon -> list of convex polygons (each an array of [x,y])
  function decompose(poly) {
    var p = cleanPoly(poly); if (p.length < 3) return null;
    if (area2(p) < 0) p.reverse();
    p = simplifyPoly(p, 40);
    if (p.length < 3) return null;
    var tris = triangulate(p), i, out = [];
    if (!tris) return [hull(p)];
    var pc = mergePieces(p, tris);
    for (i = 0; i < pc.length; i++) { var pts = []; for (var k = 0; k < pc[i].length; k++) pts.push(p[pc[i][k]]); out.push(pts); }
    return out;
  }

  // spec: {verts | parts:[convex polys] (relative to anchor), anchorWorld:[x,y], angle, isStatic}
  function makeBody(spec) {
    var polys = spec.parts ? spec.parts : [spec.verts], hs = [], infos = [], A = 0, cx = 0, cy = 0, k, i;
    for (k = 0; k < polys.length; k++) {
      var h = hull(polys[k]); if (h.length < 3) continue;
      var pr = polyProps(h); hs.push(h); infos.push(pr);
      A += pr.area; cx += pr.area * pr.cx; cy += pr.area * pr.cy;
    }
    cx /= A; cy /= A;
    var I = 0;
    for (k = 0; k < infos.length; k++) I += infos[k].inertia + infos[k].area * ((infos[k].cx - cx) * (infos[k].cx - cx) + (infos[k].cy - cy) * (infos[k].cy - cy));
    var dens = spec.density > 0 ? spec.density : 1; A *= dens; I *= dens;
    var b = { e: (spec.bounce === undefined || spec.bounce === null) ? null : spec.bounce,
              mu: (spec.friction === undefined || spec.friction === null) ? null : spec.friction, parts: [], cRel: [cx, cy], angle: spec.angle || 0, vx: 0, vy: 0, w: 0, isStatic: !!spec.isStatic,
              minx: 0, maxx: 0, miny: 0, maxy: 0 };
    b.kin = !!spec.track; b.track = spec.track || null;
    b.invM = (b.isStatic || b.kin) ? 0 : 1 / A;
    b.invI = (b.isStatic || b.kin) ? 0 : 1 / Math.max(I, 1e-6);
    for (k = 0; k < hs.length; k++) {
      var v = hs[k], n = v.length, pt = { n: n, local: [], ln: [], wv: [], wn: [], minx: 0, maxx: 0, miny: 0, maxy: 0 };
      for (i = 0; i < n; i++) pt.local.push([v[i][0] - cx, v[i][1] - cy]);
      for (i = 0; i < n; i++) {
        var j = (i + 1) % n, dx = pt.local[j][0] - pt.local[i][0], dy = pt.local[j][1] - pt.local[i][1];
        var len = Math.sqrt(dx * dx + dy * dy) || 1;
        pt.ln.push([dy / len, -dx / len]);
        pt.wv.push([0, 0]); pt.wn.push([0, 0]);
      }
      b.parts.push(pt);
    }
    var c = Math.cos(b.angle), s = Math.sin(b.angle);
    b.x = spec.anchorWorld[0] + c * cx - s * cy;
    b.y = spec.anchorWorld[1] + s * cx + c * cy;
    updateWorld(b);
    return b;
  }

  function updateWorld(b) {
    var c = Math.cos(b.angle), s = Math.sin(b.angle), i, k, p;
    b.minx = b.miny = Infinity; b.maxx = b.maxy = -Infinity;
    for (k = 0; k < b.parts.length; k++) {
      var pt = b.parts[k];
      pt.minx = pt.miny = Infinity; pt.maxx = pt.maxy = -Infinity;
      for (i = 0; i < pt.n; i++) {
        p = pt.local[i];
        var wx = b.x + c * p[0] - s * p[1], wy = b.y + s * p[0] + c * p[1];
        pt.wv[i][0] = wx; pt.wv[i][1] = wy;
        if (wx < pt.minx) pt.minx = wx; if (wx > pt.maxx) pt.maxx = wx;
        if (wy < pt.miny) pt.miny = wy; if (wy > pt.maxy) pt.maxy = wy;
        p = pt.ln[i];
        pt.wn[i][0] = c * p[0] - s * p[1]; pt.wn[i][1] = s * p[0] + c * p[1];
      }
      if (pt.minx < b.minx) b.minx = pt.minx; if (pt.maxx > b.maxx) b.maxx = pt.maxx;
      if (pt.miny < b.miny) b.miny = pt.miny; if (pt.maxy > b.maxy) b.maxy = pt.maxy;
    }
  }

  function maxSeparation(A, B) {
    var best = -Infinity, idx = 0, i, j;
    for (i = 0; i < A.n; i++) {
      var nx = A.wn[i][0], ny = A.wn[i][1], vx = A.wv[i][0], vy = A.wv[i][1], si = Infinity;
      for (j = 0; j < B.n; j++) {
        var d = nx * (B.wv[j][0] - vx) + ny * (B.wv[j][1] - vy);
        if (d < si) si = d;
      }
      if (si > best) { best = si; idx = i; }
    }
    return { sep: best, idx: idx };
  }

  function clipSeg(pts, nx, ny, off) {
    var out = [], d0 = nx * pts[0].x + ny * pts[0].y - off, d1 = nx * pts[1].x + ny * pts[1].y - off;
    if (d0 <= 0) out.push(pts[0]);
    if (d1 <= 0) out.push(pts[1]);
    if (d0 * d1 < 0) {
      var t = d0 / (d0 - d1);
      out.push({ x: pts[0].x + t * (pts[1].x - pts[0].x), y: pts[0].y + t * (pts[1].y - pts[0].y), id: 50 + (off > 0 ? 1 : 0) });
    }
    return out;
  }

  // Returns null or {nx, ny (A->B), pts:[{x,y,pen,id}]}
  function collide(A, B) {
    var sa = maxSeparation(A, B); if (sa.sep > 0) return null;
    var sb = maxSeparation(B, A); if (sb.sep > 0) return null;
    var ref, inc, ri, flip;
    if (sb.sep > sa.sep + 0.05) { ref = B; inc = A; ri = sb.idx; flip = true; }
    else { ref = A; inc = B; ri = sa.idx; flip = false; }
    var n = ref.wn[ri], r1 = ref.wv[ri], r2 = ref.wv[(ri + 1) % ref.n], i;
    var minD = Infinity, ii = 0;
    for (i = 0; i < inc.n; i++) {
      var d = n[0] * inc.wn[i][0] + n[1] * inc.wn[i][1];
      if (d < minD) { minD = d; ii = i; }
    }
    var i2 = (ii + 1) % inc.n;
    var seg = [{ x: inc.wv[ii][0], y: inc.wv[ii][1], id: ii }, { x: inc.wv[i2][0], y: inc.wv[i2][1], id: i2 }];
    var tx = r2[0] - r1[0], ty = r2[1] - r1[1], tl = Math.sqrt(tx * tx + ty * ty) || 1;
    tx /= tl; ty /= tl;
    var c1 = clipSeg(seg, -tx, -ty, -(tx * r1[0] + ty * r1[1])); if (c1.length < 2) return null;
    var c2 = clipSeg(c1, tx, ty, tx * r2[0] + ty * r2[1]); if (c2.length < 2) return null;
    var pts = [];
    for (i = 0; i < 2; i++) {
      var sep = n[0] * (c2[i].x - r1[0]) + n[1] * (c2[i].y - r1[1]);
      if (sep <= 0) pts.push({ x: c2[i].x, y: c2[i].y, pen: -sep, id: c2[i].id + (flip ? 100 : 0) });
    }
    if (!pts.length) return null;
    return { nx: flip ? -n[0] : n[0], ny: flip ? -n[1] : n[1], pts: pts };
  }

  function rectSpec(x0, y0, x1, y1) {
    return { verts: [[x0, y0], [x1, y0], [x1, y1], [x0, y1]], anchorWorld: [0, 0], angle: 0, isStatic: true };
  }

  /* cfg: { bodies:[spec], W, H, walls:{floor,left,right,ceiling}, gravity, gravityStart,
            bounce, friction, fps, frames, push, seed }
     returns array (per input body) of [{x,y,rot}] where x,y = anchor world position, rot = degrees */
  function simulate(cfg) {
    var bodies = [], i, j, k;
    for (i = 0; i < cfg.bodies.length; i++) bodies.push(makeBody(cfg.bodies[i]));
    var nUser = bodies.length, T = 4000, W = cfg.W, H = cfg.H, wl = cfg.walls || {};
    if (wl.floor)   bodies.push(makeBody(rectSpec(-T, H, W + T, H + T)));
    if (wl.ceiling) bodies.push(makeBody(rectSpec(-T, -T, W + T, 0)));
    if (wl.left)    bodies.push(makeBody(rectSpec(-T, -T, 0, H + T)));
    if (wl.right)   bodies.push(makeBody(rectSpec(W, -T, W + T, H + T)));

    for (i = 0; i < nUser; i++) {              // movable bodies without their own material use the panel defaults
      var mb = bodies[i];
      if (!mb.isStatic && !mb.kin) { if (mb.e === null) mb.e = cfg.bounce; if (mb.mu === null) mb.mu = cfg.friction; }
    }
    function centroidAt(b, x, y, deg) {
      var a = deg * Math.PI / 180, c = Math.cos(a), s = Math.sin(a);
      return [x + c * b.cRel[0] - s * b.cRel[1], y + s * b.cRel[0] + c * b.cRel[1]];
    }
    for (i = 0; i < nUser; i++) {              // velocity carried over from the layer's own keyframes
      var vf = cfg.bodies[i].velFrom, bb = bodies[i];
      if (vf && !bb.kin && !bb.isStatic && vf.dt) {
        var c1 = centroidAt(bb, vf.x, vf.y, vf.rot);
        bb.vx = (c1[0] - bb.x) / vf.dt; bb.vy = (c1[1] - bb.y) / vf.dt; bb.w = (vf.rot * Math.PI / 180 - bb.angle) / vf.dt;
      }
    }
    var seed = (cfg.seed || 1) >>> 0;
    function rnd() { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; }
    if (cfg.push) for (i = 0; i < nUser; i++) if (!bodies[i].isStatic && !bodies[i].kin) {
      bodies[i].vx += (rnd() * 2 - 1) * cfg.push; bodies[i].vy += (rnd() * 2 - 1) * cfg.push;
      bodies[i].w += (rnd() * 2 - 1) * cfg.push / 150;
    }

    var mu = cfg.friction, e = cfg.bounce, sub = Math.max(1, Math.ceil(120 / cfg.fps));
    var dt = 1 / (cfg.fps * sub), cache = {}, out = [];
    for (i = 0; i < nUser; i++) out.push([]);
    var nb = bodies.length, ITER = 10, SLOP = 0.5, BETA = 0.2, REST_V = 60;

    function record() {
      for (var q = 0; q < nUser; q++) {
        var b = bodies[q], c = Math.cos(b.angle), s = Math.sin(b.angle);
        out[q].push({ x: b.x - (c * b.cRel[0] - s * b.cRel[1]), y: b.y - (s * b.cRel[0] + c * b.cRel[1]),
                      rot: b.angle * 180 / Math.PI });
      }
    }

    function poseKinematic(b, t) {            // follow the layer's own keyframes, with matching velocity
      var tr = b.track, ft = t * cfg.fps, f0 = Math.min(Math.max(0, Math.floor(ft + 1e-9)), tr.length - 1);
      var f1 = Math.min(f0 + 1, tr.length - 1), u = f1 > f0 ? Math.min(1, ft - f0) : 0, A = tr[f0], B = tr[f1];
      var x = A.x + (B.x - A.x) * u, y = A.y + (B.y - A.y) * u, rot = A.rot + (B.rot - A.rot) * u;
      var ca = centroidAt(b, A.x, A.y, A.rot), cb = centroidAt(b, B.x, B.y, B.rot), cc = centroidAt(b, x, y, rot);
      b.x = cc[0]; b.y = cc[1]; b.angle = rot * Math.PI / 180;
      b.vx = (cb[0] - ca[0]) * cfg.fps; b.vy = (cb[1] - ca[1]) * cfg.fps; b.w = (B.rot - A.rot) * Math.PI / 180 * cfg.fps;
    }
    function step(t) {
      var g = (t >= cfg.gravityStart) ? cfg.gravity : 0, a, b, m, c;
      for (i = 0; i < nUser; i++) if (bodies[i].kin) poseKinematic(bodies[i], t);
      for (i = 0; i < nb; i++) { b = bodies[i]; if (b.invM > 0) b.vy += g * dt; updateWorld(b); }
      var cons = [], newCache = {};
      for (i = 0; i < nb; i++) for (j = i + 1; j < nb; j++) {
        a = bodies[i]; b = bodies[j];
        if (a.invM === 0 && b.invM === 0) continue;
        if (a.maxx < b.minx || b.maxx < a.minx || a.maxy < b.miny || b.maxy < a.miny) continue;
        for (var pa = 0; pa < a.parts.length; pa++) for (var pb = 0; pb < b.parts.length; pb++) {
          var P = a.parts[pa], Q = b.parts[pb];
          if (P.maxx < Q.minx || Q.maxx < P.minx || P.maxy < Q.miny || Q.maxy < P.miny) continue;
          m = collide(P, Q); if (!m) continue;
          for (k = 0; k < m.pts.length; k++) {
            var p = m.pts[k], r1x = p.x - a.x, r1y = p.y - a.y, r2x = p.x - b.x, r2y = p.y - b.y;
            var nx = m.nx, ny = m.ny, tx = ny, ty = -nx;
            var rn1 = cross2(r1x, r1y, nx, ny), rn2 = cross2(r2x, r2y, nx, ny);
            var rt1 = cross2(r1x, r1y, tx, ty), rt2 = cross2(r2x, r2y, tx, ty);
            var ee = (a.e === null) ? (b.e === null ? e : b.e) : (b.e === null ? a.e : Math.max(a.e, b.e));
            var mm = (a.mu === null) ? (b.mu === null ? mu : b.mu) : (b.mu === null ? a.mu : Math.sqrt(a.mu * b.mu));
            c = { mu: mm, a: a, b: b, nx: nx, ny: ny, tx: tx, ty: ty, r1x: r1x, r1y: r1y, r2x: r2x, r2y: r2y, pen: p.pen,
                  mN: 1 / (a.invM + b.invM + a.invI * rn1 * rn1 + b.invI * rn2 * rn2),
                  mT: 1 / (a.invM + b.invM + a.invI * rt1 * rt1 + b.invI * rt2 * rt2),
                  key: i + "_" + j + "_" + pa + "_" + pb + "_" + p.id, pn: 0, pt: 0, bias: 0 };
            var old = cache[c.key]; if (old) { c.pn = old.pn; c.pt = old.pt; }
            var dvx = b.vx - b.w * c.r2y - a.vx + a.w * c.r1y, dvy = b.vy + b.w * c.r2x - a.vy - a.w * c.r1x;
            var vn = dvx * nx + dvy * ny;
            c.bias = (vn < -REST_V) ? -ee * vn : 0;
            var Px = nx * c.pn + tx * c.pt, Py = ny * c.pn + ty * c.pt;
            a.vx -= a.invM * Px; a.vy -= a.invM * Py; a.w -= a.invI * cross2(c.r1x, c.r1y, Px, Py);
            b.vx += b.invM * Px; b.vy += b.invM * Py; b.w += b.invI * cross2(c.r2x, c.r2y, Px, Py);
            cons.push(c);
          }
        }
      }
      for (var it = 0; it < ITER; it++) for (k = 0; k < cons.length; k++) {
        c = cons[k]; a = c.a; b = c.b;
        var dx = b.vx - b.w * c.r2y - a.vx + a.w * c.r1y, dy = b.vy + b.w * c.r2x - a.vy - a.w * c.r1x;
        var dPn = c.mN * (-(dx * c.nx + dy * c.ny) + c.bias), pn0 = c.pn;
        c.pn = Math.max(pn0 + dPn, 0); dPn = c.pn - pn0;
        var Px = c.nx * dPn, Py = c.ny * dPn;
        a.vx -= a.invM * Px; a.vy -= a.invM * Py; a.w -= a.invI * cross2(c.r1x, c.r1y, Px, Py);
        b.vx += b.invM * Px; b.vy += b.invM * Py; b.w += b.invI * cross2(c.r2x, c.r2y, Px, Py);
        dx = b.vx - b.w * c.r2y - a.vx + a.w * c.r1y; dy = b.vy + b.w * c.r2x - a.vy - a.w * c.r1x;
        var dPt = c.mT * (-(dx * c.tx + dy * c.ty)), mx = c.mu * c.pn, pt0 = c.pt;
        c.pt = Math.max(-mx, Math.min(mx, pt0 + dPt)); dPt = c.pt - pt0;
        Px = c.tx * dPt; Py = c.ty * dPt;
        a.vx -= a.invM * Px; a.vy -= a.invM * Py; a.w -= a.invI * cross2(c.r1x, c.r1y, Px, Py);
        b.vx += b.invM * Px; b.vy += b.invM * Py; b.w += b.invI * cross2(c.r2x, c.r2y, Px, Py);
      }
      for (k = 0; k < cons.length; k++) newCache[cons[k].key] = { pn: cons[k].pn, pt: cons[k].pt };
      cache = newCache;
      for (i = 0; i < nb; i++) {
        b = bodies[i]; if (b.invM === 0) continue;
        var sp2 = b.vx * b.vx + b.vy * b.vy;      // safety limits for squeezed bodies
        if (sp2 > 4000 * 4000) { var kk = 4000 / Math.sqrt(sp2); b.vx *= kk; b.vy *= kk; }
        if (b.w > 25) b.w = 25; else if (b.w < -25) b.w = -25;
        b.x += b.vx * dt; b.y += b.vy * dt; b.angle += b.w * dt;
      }
      for (k = 0; k < cons.length; k++) {      // positional correction (no velocity change)
        c = cons[k]; var corr = Math.min(BETA * Math.max(c.pen - SLOP, 0), 6) / (c.a.invM + c.b.invM);
        c.a.x -= c.nx * corr * c.a.invM; c.a.y -= c.ny * corr * c.a.invM;
        c.b.x += c.nx * corr * c.b.invM; c.b.y += c.ny * corr * c.b.invM;
      }
    }

    var t = 0;
    for (var f = 0; f < cfg.frames; f++) {
      record();
      for (var s = 0; s < sub; s++) { step(t); t += dt; }
    }
    return out;
  }

  return { simulate: simulate, hull: hull, decimate: decimate, decompose: decompose };
})();
if (typeof module !== "undefined") module.exports = PhysEngine;
/* ---------------------------------------------------------------------------
   Layer glue: reads real layer geometry, runs PhysEngine, bakes keyframes.
--------------------------------------------------------------------------- */
var PhysGlue = (function () {
  var D2R = Math.PI / 180;

  function tprop(layer, mn) {
    try { return layer.property("ADBE Transform Group").property(mn); } catch (e) { return null; }
  }
  function val(prop, t, def) {
    try { var v = prop.valueAtTime(t, false); return (v === undefined || v === null) ? def : v; }
    catch (e) { try { return prop.value; } catch (e2) { return def; } }
  }
  function pv(group, mn, def) {              // value of a vector sub-property, tolerant
    try { var p = group.property(mn); return p ? p.value : def; } catch (e) { return def; }
  }
  function isStaticName(name) { return /^(static|fixed)[ _.]/i.test(name); }
  function isAnimName(name) { return /^(anim|keyed|kinematic)[ _.]/i.test(name); }


  // ---- per-layer material: Slider Controls named Bounce, Friction, Density -----------
  var MAT_NAMES = ["Bounce", "Friction", "Density"];
  function findEffect(layer, name) {
    try {
      var fx = layer.property("ADBE Effect Parade");
      for (var i = 1; i <= fx.numProperties; i++) {
        var f = fx.property(i);
        if (String(f.name).toLowerCase() === name.toLowerCase()) return f;
      }
    } catch (e) {}
    return null;
  }
  function readMaterial(layer, t0) {
    var m = {}, i;
    for (i = 0; i < MAT_NAMES.length; i++) {
      var f = findEffect(layer, MAT_NAMES[i]);
      if (!f) continue;
      try { var v = f.property(1).valueAtTime(t0, false); if (typeof v === "number" && isFinite(v)) m[MAT_NAMES[i].toLowerCase()] = v; } catch (e) {}
    }
    return m;
  }
  function addControlsTo(layer, defaults) {      // returns number of controls added
    var fx = layer.property("ADBE Effect Parade"), added = 0, i;
    for (i = 0; i < MAT_NAMES.length; i++) {
      if (findEffect(layer, MAT_NAMES[i])) continue;
      var f = fx.addProperty("ADBE Slider Control");
      f.name = MAT_NAMES[i];
      f.property(1).setValue(defaults[MAT_NAMES[i].toLowerCase()]);
      added++;
    }
    return added;
  }
  function addControls(defaults) {
    var comp = app.project.activeItem;
    if (!comp || !(comp instanceof CompItem)) throw new Error("Open a composition first.");
    var sel = comp.selectedLayers; if (!sel || !sel.length) throw new Error("Select the layers first.");
    var n = 0, i;
    app.beginUndoGroup("Add Physics Controls");
    try { for (i = 0; i < sel.length; i++) n += addControlsTo(sel[i], defaults); } finally { app.endUndoGroup(); }
    return "Added " + n + " control(s) (Bounce, Friction, Density) to " + sel.length + " layer(s). Edit them in each layer's Effects. Select several layers to edit them together.";
  }

  // ---- shape layer: gather points in layer space -----------------------------------
  function groupXf(tr, parent) {
    var pos = pv(tr, "ADBE Vector Position", [0, 0]), anc = pv(tr, "ADBE Vector Anchor", [0, 0]);
    var sc = pv(tr, "ADBE Vector Scale", [100, 100]), rot = pv(tr, "ADBE Vector Rotation", 0) * D2R;
    var c = Math.cos(rot), s = Math.sin(rot);
    return function (p) {
      var x = (p[0] - anc[0]) * sc[0] / 100, y = (p[1] - anc[1]) * sc[1] / 100;
      return parent([c * x - s * y + pos[0], s * x + c * y + pos[1]]);
    };
  }
  function addBezier(sh, xf, out) {
    var v = sh.vertices, it = sh.inTangents || [], ot = sh.outTangents || [], n = v.length, i, k;
    for (i = 0; i < n; i++) {
      var j = (i + 1) % n; if (j === 0 && sh.closed === false) { out.push(xf(v[i])); break; }
      var o = ot[i] || [0, 0], inn = it[j] || [0, 0];
      var straight = !o[0] && !o[1] && !inn[0] && !inn[1];
      out.push(xf(v[i]));
      if (!straight) for (k = 1; k < 10; k++) {
        var u = k / 10, w = 1 - u;
        var p1 = [v[i][0] + o[0], v[i][1] + o[1]], p2 = [v[j][0] + inn[0], v[j][1] + inn[1]];
        out.push(xf([w * w * w * v[i][0] + 3 * w * w * u * p1[0] + 3 * w * u * u * p2[0] + u * u * u * v[j][0],
                     w * w * w * v[i][1] + 3 * w * w * u * p1[1] + 3 * w * u * u * p2[1] + u * u * u * v[j][1]]));
      }
    }
  }
  // each shape item becomes its own closed polygon (array of [x,y]) in layer space
  function collectShape(group, xf, polys) {
    var n = group.numProperties, i, pp;
    for (i = 1; i <= n; i++) {
      var p = group.property(i), mn = p.matchName;
      if (mn === "ADBE Vector Group") {
        collectShape(p.property("ADBE Vectors Group"), groupXf(p.property("ADBE Vector Transform Group"), xf), polys);
      } else if (mn === "ADBE Vector Shape - Group") {
        pp = []; addBezier(p.property("ADBE Vector Shape").value, xf, pp); polys.push(pp);
      } else if (mn === "ADBE Vector Shape - Rect") {
        var rs = pv(p, "ADBE Vector Rect Size", [100, 100]), rp = pv(p, "ADBE Vector Rect Position", [0, 0]);
        polys.push([xf([rp[0] - rs[0] / 2, rp[1] - rs[1] / 2]), xf([rp[0] + rs[0] / 2, rp[1] - rs[1] / 2]),
                    xf([rp[0] + rs[0] / 2, rp[1] + rs[1] / 2]), xf([rp[0] - rs[0] / 2, rp[1] + rs[1] / 2])]);
      } else if (mn === "ADBE Vector Shape - Ellipse") {
        var es = pv(p, "ADBE Vector Ellipse Size", [100, 100]), ep = pv(p, "ADBE Vector Ellipse Position", [0, 0]);
        pp = [];
        for (var k = 0; k < 24; k++) { var a = k * Math.PI / 12; pp.push(xf([ep[0] + Math.cos(a) * es[0] / 2, ep[1] + Math.sin(a) * es[1] / 2])); }
        polys.push(pp);
      } else if (mn === "ADBE Vector Shape - Star") {
        var cnt = Math.max(3, Math.round(pv(p, "ADBE Vector Star Points", 5))), orad = pv(p, "ADBE Vector Star Outer Radius", 50);
        var irad = pv(p, "ADBE Vector Star Inner Radius", orad * 0.5), isStar = pv(p, "ADBE Vector Star Type", 1) === 1;
        var sp = pv(p, "ADBE Vector Star Position", [0, 0]), srot = pv(p, "ADBE Vector Star Rotation", 0) * D2R;
        pp = [];
        for (var q = 0; q < cnt; q++) {
          var aa = srot + q * 2 * Math.PI / cnt - Math.PI / 2;
          pp.push(xf([sp[0] + Math.cos(aa) * orad, sp[1] + Math.sin(aa) * orad]));
          if (isStar) { var ab = aa + Math.PI / cnt; pp.push(xf([sp[0] + Math.cos(ab) * irad, sp[1] + Math.sin(ab) * irad])); }
        }
        polys.push(pp);
      }
    }
  }
  function shapePolys(layer) {
    var root = null, polys = [], total = 0, i;
    try { root = layer.property("ADBE Root Vectors Group"); } catch (e) {}
    if (!root) return null;
    try { collectShape(root, function (p) { return p; }, polys); } catch (e) { return null; }
    var good = [];
    for (i = 0; i < polys.length; i++) if (polys[i].length >= 3) { good.push(polys[i]); total += polys[i].length; }
    return good.length ? good : null;
  }

  // ---- build a collider spec for one layer ---------------------------------------
  // mode: 0 = follow outline (concave OK), 1 = convex hull, 2 = bounding box
  function buildSpec(layer, t0, mode) {
    var warn = null;
    try { if (layer.parent) return { skip: "has a parent (unparent it first)" }; } catch (e) {}
    try { if (layer.threeDLayer) return { skip: "is 3D (2D layers only)" }; } catch (e) {}
    var pos = val(tprop(layer, "ADBE Position"), t0, [0, 0]);
    var anc = val(tprop(layer, "ADBE Anchor Point"), t0, [0, 0]);
    var sc = val(tprop(layer, "ADBE Scale"), t0, [100, 100]);
    var rot = val(tprop(layer, "ADBE Rotate Z"), t0, 0);
    function toRel(pt) { return [(pt[0] - anc[0]) * sc[0] / 100, (pt[1] - anc[1]) * sc[1] / 100]; }
    var polys = (mode !== 2) ? shapePolys(layer) : null, spec = null, how = "bounding box", i, k;
    if (polys) {
      var rels = [];
      for (i = 0; i < polys.length; i++) { var rp = []; for (k = 0; k < polys[i].length; k++) rp.push(toRel(polys[i][k])); rels.push(rp); }
      if (mode === 0) {
        var parts = [];
        for (i = 0; i < rels.length; i++) { var d = PhysEngine.decompose(rels[i]); if (d) for (k = 0; k < d.length; k++) if (d[k].length >= 3) parts.push(d[k]); }
        if (parts.length) { spec = { parts: parts }; how = "shape outline, " + parts.length + " convex piece" + (parts.length > 1 ? "s" : ""); }
      } else {
        var all = []; for (i = 0; i < rels.length; i++) all = all.concat(rels[i]);
        var h = PhysEngine.decimate(PhysEngine.hull(all), 24);
        if (h.length >= 3) { spec = { verts: h }; how = "convex hull"; }
      }
    }
    if (!spec) {
      var r = null;
      try { r = layer.sourceRectAtTime(t0, false); } catch (e) {}
      if (!r || !(r.width > 0) || !(r.height > 0)) {
        var w = 0, hh = 0;
        try { w = layer.width; hh = layer.height; } catch (e2) {}
        if (!(w > 0) || !(hh > 0)) return { skip: "has no measurable bounds" };
        r = { left: 0, top: 0, width: w, height: hh };
        warn = "used layer width/height";
      }
      var bp = [[r.left, r.top], [r.left + r.width, r.top], [r.left + r.width, r.top + r.height], [r.left, r.top + r.height]];
      var br = []; for (i = 0; i < 4; i++) br.push(toRel(bp[i]));
      spec = { verts: br };
    }
    spec.anchorWorld = [pos[0], pos[1]]; spec.angle = rot * D2R; spec.isStatic = isStaticName(layer.name);
    var mat = readMaterial(layer, t0);
    if (mat.bounce !== undefined) spec.bounce = Math.max(0, Math.min(1.5, mat.bounce));
    if (mat.friction !== undefined) spec.friction = Math.max(0, mat.friction);
    if (mat.density !== undefined && mat.density > 0) spec.density = mat.density;
    return { spec: spec, posZ: pos.length > 2 ? pos[2] : null, how: how, warn: warn };
  }

  function removeKeysFrom(prop, t0) {
    if (!prop) return;
    try { prop.expressionEnabled = false; } catch (e) {}
    for (var i = prop.numKeys; i >= 1; i--) if (prop.keyTime(i) >= t0 - 1e-6) prop.removeKey(i);
  }


  // ---- split layers into independent pieces ----------------------------------------
  function isTextLayer(L) { try { return !!L.property("ADBE Text Properties"); } catch (e) { return false; } }
  function setStatic(prop, v) {
    try { prop.expressionEnabled = false; } catch (e) {}
    while (prop.numKeys > 0) prop.removeKey(1);
    prop.setValue(v);
  }
  function splitShapeGroups(L, log) {
    var root = L.property("ADBE Root Vectors Group"), gs = [], i, made = [];
    for (i = 1; i <= root.numProperties; i++) if (root.property(i).matchName === "ADBE Vector Group") gs.push(i);
    if (gs.length < 2) { log("'" + L.name + "': only one shape in the layer, nothing to split."); return made; }
    for (var k = 0; k < gs.length; k++) {
      var dup = L.duplicate(), r2 = dup.property("ADBE Root Vectors Group"), list = [];
      for (i = 1; i <= r2.numProperties; i++) if (r2.property(i).matchName === "ADBE Vector Group") list.push(r2.property(i));
      for (i = list.length - 1; i >= 0; i--) if (i !== k) list[i].remove();
      dup.name = L.name + "_" + (k + 1);
      made.push(dup);
    }
    L.enabled = false;
    log("'" + L.name + "': split into " + made.length + " shape layers.");
    return made;
  }
  function setText(tmp, str, left) {
    var d = tmp.property("ADBE Text Properties").property("ADBE Text Document"), v = d.value;
    v.text = str;
    if (left) { try { v.justification = ParagraphJustification.LEFT_JUSTIFY; } catch (e) {} }
    d.setValue(v);
  }
  function splitText(L, t0, log) {
    var doc = L.property("ADBE Text Properties").property("ADBE Text Document"), td = doc.value, made = [];
    var lines = String(td.text).split(/\r\n|\r|\n/), just = null;
    try { just = td.justification; } catch (e) {}
    var isCenter = false, isRight = false;
    try { isCenter = (just === ParagraphJustification.CENTER_JUSTIFY); isRight = (just === ParagraphJustification.RIGHT_JUSTIFY); } catch (e) {}
    var tmp = L.duplicate(); tmp.name = "__physics_measure__";
    function right(str) { setText(tmp, str + "|", true); var r = tmp.sourceRectAtTime(t0, false); return r.left + r.width; }
    function bottom(str) { setText(tmp, str, true); var r = tmp.sourceRectAtTime(t0, false); return r.top + r.height; }
    var items = [];
    try {
      var bar = right("");
      if (!(bar > 0)) throw new Error("This build cannot measure text. Convert the text to shapes first (Layer > Create Shapes from Text), then press Split again.");
      var barBottom = bottom("|"), j, i;
      for (j = 0; j < lines.length; j++) {
        var line = lines[j], adv = right(line) - bar;
        var startX = isCenter ? -adv / 2 : (isRight ? -adv : 0);
        var yOff = 0;
        if (j > 0) yOff = bottom(lines.slice(0, j).join("\r") + "\r|") - barBottom;
        for (i = 0; i < line.length; i++) {
          var ch = line.charAt(i); if (/\s/.test(ch)) continue;
          items.push({ ch: ch, dx: startX + (i > 0 ? right(line.substring(0, i)) - bar : 0), dy: yOff });
        }
      }
    } finally { tmp.remove(); }
    if (!items.length) { log("'" + L.name + "': no visible letters."); return made; }
    var pos = val(tprop(L, "ADBE Position"), t0, [0, 0]), sc = val(tprop(L, "ADBE Scale"), t0, [100, 100]);
    var rot = val(tprop(L, "ADBE Rotate Z"), t0, 0) * D2R, c = Math.cos(rot), s = Math.sin(rot);
    for (var n = 0; n < items.length; n++) {
      var it = items[n], dup = L.duplicate(), lx = it.dx * sc[0] / 100, ly = it.dy * sc[1] / 100;
      setText(dup, it.ch, true);
      var np = [pos[0] + c * lx - s * ly, pos[1] + s * lx + c * ly]; if (pos.length > 2) np.push(pos[2]);
      setStatic(tprop(dup, "ADBE Position"), np);
      dup.name = L.name + "_" + (n + 1);
      made.push(dup);
    }
    L.enabled = false;
    log("'" + L.name + "': split into " + made.length + " letter layers.");
    return made;
  }
  function split(defaults, log) {
    var comp = app.project.activeItem;
    if (!comp || !(comp instanceof CompItem)) throw new Error("Open a composition first.");
    var src = [], sel = comp.selectedLayers, i, all = [];
    if (!sel || !sel.length) throw new Error("Select a text layer, or a shape layer that contains several shapes.");
    for (i = 0; i < sel.length; i++) src.push(sel[i]);
    app.beginUndoGroup("Split Into Pieces");
    try {
      for (i = 0; i < src.length; i++) {
        var made = [];
        try {
          if (isTextLayer(src[i])) made = splitText(src[i], comp.time, log);
          else made = splitShapeGroups(src[i], log);
        } catch (err) { log("'" + src[i].name + "': " + (err && err.message ? err.message : err)); }
        if (made.length) { src[i].selected = false; for (var k = 0; k < made.length; k++) { made[k].selected = true; all.push(made[k]); } }
      }
      for (i = 0; i < all.length; i++) { try { addControlsTo(all[i], defaults); } catch (e) {} }
    } finally { app.endUndoGroup(); }
    return all.length ? all.length + " new layers created and selected. Press Simulate." : "Nothing was split.";
  }

  function sampleTrack(layer, t0, fps, frames) {
    var pp = tprop(layer, "ADBE Position"), rp = tprop(layer, "ADBE Rotate Z"), tr = [], f;
    for (f = 0; f < frames; f++) {
      var t = t0 + f / fps, p = val(pp, t, [0, 0]);
      tr.push({ x: p[0], y: p[1], rot: val(rp, t, 0) });
    }
    return tr;
  }
  function isAnimated(layer) {
    try { var pp = tprop(layer, "ADBE Position"), rp = tprop(layer, "ADBE Rotate Z");
      return (pp && (pp.numKeys > 0 || pp.expressionEnabled)) || (rp && (rp.numKeys > 0 || rp.expressionEnabled)); } catch (e) { return false; }
  }

  // opts: {gravityOn, gravity, gravityStart, bounce, friction, duration, push, carry, floor,left,right,ceiling, mode}
  function bake(opts, log) {
    var comp = app.project.activeItem;
    if (!comp || !(comp instanceof CompItem)) throw new Error("Open a composition first.");
    var sel = comp.selectedLayers, selIdx = {}, i;
    if (!sel || !sel.length) throw new Error("Select the layers you want to simulate.");
    for (i = 0; i < sel.length; i++) selIdx[sel[i].index] = true;
    var t0 = comp.time, fps = comp.frameRate, specs = [], layers = [];
    var frames = Math.max(2, Math.round(opts.duration * fps) + 1);
    if (frames > 20000) throw new Error("Duration too long.");
    for (i = 1; i <= comp.numLayers; i++) {
      var L = comp.layer(i), isSt = isStaticName(L.name), isKin = isAnimName(L.name);
      if (!selIdx[i] && !isSt && !isKin) continue;
      var r = buildSpec(L, t0, opts.mode);
      if (r.skip) { log("Skipped '" + L.name + "': " + r.skip); continue; }
      var tag = "";
      if (isKin) { r.spec.track = sampleTrack(L, t0, fps, frames); tag = "[animated obstacle] "; }
      else if (isSt) tag = "[static] ";
      else if (opts.carry && isAnimated(L)) {
        var h = 1 / fps, t1 = (t0 >= h - 1e-9) ? t0 - h : t0 + h, p1 = val(tprop(L, "ADBE Position"), t1, r.spec.anchorWorld);
        r.spec.velFrom = { x: p1[0], y: p1[1], rot: val(tprop(L, "ADBE Rotate Z"), t1, r.spec.angle / D2R), dt: t1 - t0 };
        tag = "[keeps its keyframed motion] ";
      }
      specs.push(r.spec); layers.push({ layer: L, info: r });
      log(tag + L.name + ": " + r.how + (r.warn ? " (" + r.warn + ")" : ""));
    }
    var dyn = 0; for (i = 0; i < specs.length; i++) if (!specs[i].isStatic && !specs[i].track) dyn++;
    if (!dyn) throw new Error("No movable layers selected (layers named Static/Fixed/Anim are obstacles).");
    var res = PhysEngine.simulate({
      bodies: specs, W: comp.width, H: comp.height,
      walls: { floor: opts.floor, left: opts.left, right: opts.right, ceiling: opts.ceiling },
      gravity: opts.gravityOn ? opts.gravity : 0, gravityStart: opts.gravityStart,
      bounce: opts.bounce, friction: opts.friction, fps: fps, frames: frames, push: opts.push, seed: 7
    });
    app.beginUndoGroup("Bake 2D Physics");
    try {
      for (i = 0; i < layers.length; i++) {
        if (specs[i].isStatic || specs[i].track) continue;
        var L2 = layers[i].layer, pp = tprop(L2, "ADBE Position"), rp = tprop(L2, "ADBE Rotate Z"), z = layers[i].info.posZ;
        removeKeysFrom(pp, t0); removeKeysFrom(rp, t0);
        for (var f = 0; f < frames; f++) {
          var s = res[i][f], t = t0 + f / fps;
          pp.setValueAtTime(t, z === null ? [s.x, s.y] : [s.x, s.y, z]);
          rp.setValueAtTime(t, s.rot);
        }
      }
    } finally { app.endUndoGroup(); }
    return dyn + " layer(s) baked, " + frames + " frames from " + t0.toFixed(2) + "s.";
  }

  function clear() {
    var comp = app.project.activeItem;
    if (!comp || !(comp instanceof CompItem)) throw new Error("Open a composition first.");
    var sel = comp.selectedLayers, t0 = comp.time;
    app.beginUndoGroup("Clear Physics Keys");
    try {
      for (var i = 0; i < sel.length; i++) { removeKeysFrom(tprop(sel[i], "ADBE Position"), t0); removeKeysFrom(tprop(sel[i], "ADBE Rotate Z"), t0); }
    } finally { app.endUndoGroup(); }
    return "Cleared position/rotation keys from " + t0.toFixed(2) + "s on " + sel.length + " layer(s).";
  }
  return { bake: bake, clear: clear, buildSpec: buildSpec, split: split, addControls: addControls };
})();
if (typeof module !== "undefined") module.exports = PhysGlue;
/* ---------------------------------------------------------------------------
   Panel UI (ScriptUI). Runs as a dockable panel or a floating palette.
--------------------------------------------------------------------------- */
(function (thisObj) {
  function num(et, def) { var v = parseFloat(et.text); return isNaN(v) ? def : v; }
  function row(parent, label, value) {
    var g = parent.add("group"); g.orientation = "row";
    var s = g.add("statictext", undefined, label); s.preferredSize = [130, 20];
    var e = g.add("edittext", undefined, value); e.preferredSize = [70, 22];
    return e;
  }
  var isPanel = (typeof Panel !== "undefined") && (thisObj instanceof Panel);
  var win = isPanel ? thisObj : new Window("palette", "2D Physics", undefined, { resizeable: true });
  win.orientation = "column"; win.alignChildren = ["fill", "top"]; win.spacing = 6; win.margins = 10;

  var cGrav = win.add("checkbox", undefined, "Gravity enabled"); cGrav.value = true;
  var cCarry = win.add("checkbox", undefined, "Carry velocity from keyframes"); cCarry.value = true;
  var eGrav = row(win, "Gravity (px/s\u00B2)", "1500");
  var eStart = row(win, "Gravity starts at (s)", "0");
  var eBounce = row(win, "Default bounce (0-1)", "0.35");
  var eFric = row(win, "Default friction (0-1)", "0.5");
  var ePush = row(win, "Initial push (px/s)", "0");
  var eDur = row(win, "Duration (s)", "6");

  var wg = win.add("group"); wg.orientation = "row";
  wg.add("statictext", undefined, "Walls:");
  var cF = wg.add("checkbox", undefined, "Floor"); cF.value = true;
  var cL = wg.add("checkbox", undefined, "Left"); cL.value = true;
  var cR = wg.add("checkbox", undefined, "Right"); cR.value = true;
  var cC = wg.add("checkbox", undefined, "Top"); cC.value = false;

  var mg = win.add("group"); mg.orientation = "row";
  mg.add("statictext", undefined, "Collider:");
  var dd = mg.add("dropdownlist", undefined, ["Follow shape outline (concave)", "Convex hull", "Bounding box"]); dd.selection = 0;

  var bGo = win.add("button", undefined, "Simulate selected layers");
  var bSplit = win.add("button", undefined, "Split selected into letters / shapes");
  var bCtl = win.add("button", undefined, "Add Bounce/Friction/Density controls");
  var bClr = win.add("button", undefined, "Clear keys from current time");
  var st = win.add("statictext", undefined, "Select layers, put the playhead where the sim should start.", { multiline: true });
  st.preferredSize = [280, 90];
  win.add("statictext", undefined, "Layer names (the prefix must be followed by a space, _ or .):\n  Static_ / Fixed_ = fixed obstacle\n  Anim_ / Keyed_ / Kinematic_ = follows its own keyframes\nAnything else selected = moved by physics.\nPer-layer Bounce / Friction / Density: use the controls button, then edit each layer's sliders.", { multiline: true }).preferredSize = [280, 110];

  var lines = [];
  function log(m) { lines.push(m); }
  function run(fn) {
    lines = [];
    try { var msg = fn(); lines.push(msg); } catch (err) { lines.push("Error: " + (err && err.message ? err.message : err)); }
    st.text = lines.join("\n");
  }
  bGo.onClick = function () {
    run(function () {
      return PhysGlue.bake({
        gravityOn: cGrav.value, gravity: num(eGrav, 1500), gravityStart: Math.max(0, num(eStart, 0)),
        bounce: Math.max(0, Math.min(1, num(eBounce, 0.35))), friction: Math.max(0, num(eFric, 0.5)),
        push: Math.max(0, num(ePush, 0)), carry: cCarry.value, duration: Math.max(0.1, num(eDur, 6)),
        floor: cF.value, left: cL.value, right: cR.value, ceiling: cC.value,
        mode: dd.selection ? dd.selection.index : 0
      }, log);
    });
  };
  bClr.onClick = function () { run(function () { return PhysGlue.clear(); }); };
  function defaults() { return { bounce: Math.max(0, Math.min(1, num(eBounce, 0.35))), friction: Math.max(0, num(eFric, 0.5)), density: 1 }; }
  bSplit.onClick = function () { run(function () { return PhysGlue.split(defaults(), log); }); };
  bCtl.onClick = function () { run(function () { return PhysGlue.addControls(defaults()); }); };

  if (win.layout) { win.layout.layout(true); }
  if (!isPanel) { try { win.center(); } catch (e) {} win.show(); }
})(this);
