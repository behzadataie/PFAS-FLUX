# Monte Carlo design

The uncertainty calculation is a transparent screening calculation for the technical note benchmark. It is not a formal Bayesian calibration or a replacement for site-specific uncertainty analysis. Each realization samples uncertain input quantities, recomputes the same mass-discharge equations used in the deterministic run, and stores percentile intervals and exceedance probabilities.

## Sampled parameters and distributions

| Parameter | Range source | Distribution | Reason |
| --- | --- | --- | --- |
| Hydraulic conductivity K_i | assumed_K_zones.csv: K_low_m_d, K_central_m_d, K_high_m_d | Log-triangular | K is positive and commonly varies by orders of magnitude. Sampling in log10 space prevents the high end of K from dominating only because arithmetic units were chosen. |
| Hydraulic gradient I_i | assumed_K_zones.csv plus control_plane_cells.csv: gradient_normal_central | Triangular | Low and high ranges come from K-zone assumptions; the central value uses the control-plane normal gradient when supplied. |
| Effective saturated thickness | assumed_K_zones.csv: effective_thickness_low_m, central, high | Triangular | The Monte Carlo recomputes cell area as cell width times sampled thickness. |
| Concentration C_i,k | control_plane_cells.csv or surface_water_nodes.csv or wells.csv central concentration; CLI multipliers default to 0.5 and 2.0 | Triangular | Low = central x concentration_low_multiplier; mode = central; high = central x concentration_high_multiplier. Zero or missing values remain zero. |
| Completeness factor F_comp,G | analyte_groups.csv: F_comp_low, F_comp_central, F_comp_high | Triangular | Scenario factor for analyte incompleteness or precursor/TOP uncertainty; not a measurement. |
| Surface-water/event flow | surface_water_nodes.csv: flow_low_m3_d, flow_central_m3_d, flow_high_m3_d | Triangular | Produces g/d if the flow input is m3/d, or g/event if the user supplies an event volume on the same numeric basis. |
| Receptor mixing flow Q_mix | receptor_flows.csv: Q_low_m3_d, Q_central_m3_d, Q_high_m3_d | Triangular | Represents receiving-water flow, extraction flow, or decision-specific mixing flow. |
| Mixing/allocation factor f_mix | receptor_flows.csv: mixing_factor_low, central, high | Triangular | Represents the fraction of receptor flow assigned to the source contribution. |

## Calculation sequence in each realization

1. Sample K_i, hydraulic gradient, and effective saturated thickness for each control-plane cell.
2. Recompute cell area where uncertain thickness is sampled.
3. Sample each central concentration with the default triangular multiplier range of 0.5 to 2.0, unless changed by command-line options.
4. Compute q_i = K_i I_i.
5. Compute M_d,k = sum_i(C_i,k q_i A_i) x 1e-3.
6. Sum included analytes to obtain measured group load M_d,G.
7. Sample F_comp,G and compute M*_d,G = F_comp,G M_d,G.
8. Sample receptor flow Q_mix and mixing factor f_mix. Convert the placeholder criterion from ng/L to ug/L.
9. Compute M_allow = C_crit Q_mix f_mix x 1e-3.
10. Compute P_exceed as the fraction of realizations where current load exceeds allowable load.

## Random seeds

The default seed is 20260703. Separate random-number streams are created by adding deterministic offsets for surface-water, well and receptor calculations. This makes the run reproducible while keeping the sampling streams distinct.

## Independence assumption

The default implementation samples input variables independently. This is deliberate for a short technical note benchmark, but real sites may require correlated K fields, concentration interpolation uncertainty, time-series uncertainty, and censored-data treatment.

## Defaults used in the provided run

```text
selected_group_id = G3
n_realizations = 10000
seed = 20260703
concentration_low_multiplier = 0.5
concentration_high_multiplier = 2.0
receptor_load_type = F_comp adjusted scenario
```
