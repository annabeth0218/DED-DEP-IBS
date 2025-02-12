o1 <- magma.deddep |>
  filter(P_ded < .02) |>
  filter(P_dep < .02)

nrow(o1)

o2 <- magma.ibsdep |>
  filter(P_ibs < .01) |>
  filter(P_dep < .01)

nrow(o2)

o3 <- magma.ibsded |>
  filter(P_ibs < .01) |>
  filter(P_ded < .01)

nrow(o3)

magma.ov3 <- inner_join(magma.ibsded, magma_dep, by = "GENE") |>
  select(GENE, SYMBOL = SYMBOL.x, CHR = CHR.x, START = START.x,
         STOP = STOP.x, Nsnps_ibs = Nsnps_ibs,
         Nsnps_ded = Nsnps_ded, Nsnps_dep = NSNPS,
         P_ibs = P_ibs, P_ded = P_ded, P_dep = P)


