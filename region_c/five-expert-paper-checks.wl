(* Exact moving-endpoint checks for the final two Region C traces.
   No Integrate, NIntegrate, or hypergeometric simplification is used.
   This checks identities, not the analytic endpoint estimates in the paper. *)
ClearAll[q, s, x, sp, cc, pp, ff, clean, step, rep, jj, kk, ii];
$Assumptions = 0 < q < 1 && 0 < x < 1 && q <= s <= 1;
clean[e_List] := clean /@ e;
clean[e_] := Factor[Together[e]];
sp[x_] := x^6 - 3 x^5 - 3 x^4 + 12 x^3 Log[x] + 3 x^2 + 3 x - 1;
cc[q_] := -2 Sqrt[2] q^(3/2) sp[q]/(9 (q - 1)^8 (q + 1));
pp[0] = 4 q^4 x + 2 q^4 + 3 q^3 x^2 + 2 q^3 x + q^3
 + s (-6 q^3 x^2 - 4 q^3 x - 2 q^3 - 2 q x^4 - 4 q x^3 - 6 q x^2)
 + s^2 (q x^4 + 2 q x^3 + 3 q x^2 + 2 x^4 + 4 x^3);
pp[1] = 10 q^4 x + 5 q^4 + 3 q^3 x^2 + 2 q^3 x + q^3
 + 9 q^2 (x^3 + x^2 + x)
 + s (3 q^3 x^2 + 2 q^3 x + q^3 - 36 q^2 (x^3 + x^2 + x)
      + q x^4 + 2 q x^3 + 3 q x^2)
 + s^2 (9 q^2 (x^3 + x^2 + x) + q x^4 + 2 q x^3 + 3 q x^2
        + 5 x^4 + 10 x^3);
ff[0] = -Sqrt[2] sp[x] pp[0]/(9 q^2 Sqrt[x] (x - 1)^9 (x + 1)^2);
ff[1] = -Sqrt[2] sp[x] pp[1]/(36 q^2 Sqrt[x] (x - 1)^9 (x + 1)^2);

(* {a,b} denotes a(q,s) + Integral[b(q,s,x),{x,0,q}].
   In D[b,q], x is independent; only afterward set x->q. *)
rep[j_] := {cc[q] (s - 1)^2, ff[j]};
step[{a_, b_}] := {
 clean[D[a, q] + (b /. x -> q)], clean[D[b, q]]};
jj[j_] := Map[clean[-q^2 D[#, {s, 2}]] &, rep[j]];
kk[j_] := Map[clean[q^2 (# /. s -> 1)] &, rep[j]];
ii[j_] := Map[clean[# /. s -> q] &, rep[j]];

b0 = 6 q^5 - 60 q^4 Log[q] + 125 q^4 - 120 q^3 Log[q]
 - 80 q^3 - 60 q^2 + 10 q - 1;
b1 = 22 q^6 - 300 q^5 Log[q] + 763 q^5 - 1260 q^4 Log[q]
 + 615 q^4 - 840 q^3 Log[q] - 1180 q^3 - 250 q^2 + 33 q - 3;
m0 = q^8 - 11 q^7 + 104 q^6 - 540 q^5 Log[q] + 1021 q^5
 - 1320 q^4 Log[q] - 540 q^3 Log[q] - 1021 q^3 - 104 q^2 + 11 q - 1;
m1 = q^8 - 11 q^7 + 184 q^6 - 1500 q^5 Log[q] + 3261 q^5
 - 4200 q^4 Log[q] - 1500 q^3 Log[q] - 3261 q^3 - 184 q^2 + 11 q - 1;
n0 = q^6 - 24 q^5 + 180 q^4 Log[q] - 375 q^4 + 480 q^3 Log[q]
 + 180 q^2 Log[q] + 375 q^2 + 24 q - 1;
n1 = q^6 - 34 q^5 + 300 q^4 Log[q] - 655 q^4 + 840 q^3 Log[q]
 + 300 q^2 Log[q] + 655 q^2 + 34 q - 1;

expected = {
 -2 Sqrt[2] q^(3/2) b0/(3 (q - 1)^10),
 Sqrt[2] Sqrt[q] b1/(2 (q - 1)^11),
 -3 Sqrt[2] m0/(4 q^(3/2) (q - 1)^11),
 -Sqrt[2] m1/(4 q^(3/2) (q - 1)^11),
 Sqrt[2] n0/(2 Sqrt[q] (q - 1)^9),
 Sqrt[2] n1/(4 Sqrt[q] (q - 1)^9)};
computed = {Nest[step, jj[0], 2], Nest[step, jj[1], 3],
 Nest[step, kk[0], 5], Nest[step, kk[1], 5],
 Nest[step, ii[0], 3], Nest[step, ii[1], 3]};
identityChecks = MapThread[
  Function[{pair, target}, clean[pair - {target, 0}] === {0, 0}],
  {computed, expected}];
Print["J0'', J1''', K0^(5), K1^(5), I0''', I1''': ", identityChecks];

signChecks = {
 clean[D[b0, {q, 5}] - 720 (q - 1)^2/q^2] === 0,
 clean[D[b1, {q, 6}] - 1440 (q - 1) (11 q^2 - 14 q + 7)/q^3] === 0,
 clean[D[m0/(60 q^3 (9 q^2 + 22 q + 9)), q]
   - (q - 1)^6 (27 q^4 + 52 q^3 + 162 q^2 + 52 q + 27)
     /(60 q^4 (9 q^2 + 22 q + 9)^2)] === 0,
 clean[D[m1/(300 q^3 (5 q^2 + 14 q + 5)), q]
   - (q - 1)^6 (5 q^4 + 12 q^3 + 158 q^2 + 12 q + 5)
     /(100 q^4 (5 q^2 + 14 q + 5)^2)] === 0,
 clean[D[n0/(60 q^2 (3 q^2 + 8 q + 3)), q]
   - (q - 1)^8/(10 q^3 (3 q^2 + 8 q + 3)^2)] === 0,
 clean[D[n1/(60 q^2 (5 q^2 + 14 q + 5)), q]
   - (q - 1)^6 (5 q^2 - 34 q + 5)
     /(30 q^3 (5 q^2 + 14 q + 5)^2)] === 0};
Print["Terminal sign identities: ", signChecks];
If[!(And @@ Join[identityChecks, signChecks]), Print["CHECK FAILED"]; Quit[1]];
Print["PASS: all twelve exact derivative identities."];
