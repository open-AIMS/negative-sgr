## Arm colours for every figure in the study, validated rather than chosen by eye.
##
## Eight series is past the point where a categorical palette is safe by
## inspection, so these were run through the dataviz validator in both modes
## (lightness band, chroma floor, colour-vision separation, normal-vision floor,
## contrast against the surface).
##
## Two things the validation forced:
##
## 1. `arm_levels` is the LEGEND ORDER, and it is not arbitrary. The obvious
##    grouping -- retained, imposed, family -- puts D (violet) next to B1 (blue),
##    which fails deutan separation at delta-E 1.9. Reversing the middle group so
##    the boundary falls at D-B2 clears every adjacent pair in both modes.
## 2. All-pairs separation is NOT achievable with eight hues; adjacent-pair is.
##    That is safe only because the figures facet by `arm_group`, so series that
##    are close in hue never share a panel. Any figure that puts all eight in one
##    panel must re-validate with `--pairs all` first, and will fail.
arm_levels <- c("A", "C", "D", "B2", "B3", "B1", "E", "F")

arm_group <- c(A = "Measurement retained", C = "Measurement retained",
               D = "Measurement retained",
               B1 = "Zero boundary imposed", B2 = "Zero boundary imposed",
               B3 = "Zero boundary imposed",
               E = "Zero bounded by the family",
               F = "Zero bounded by the family")

arm_cols_light <- c(A = "#0E7C61", C = "#C25E12", D = "#8E44C9",
                    B2 = "#7C7A00", B3 = "#B0217A", B1 = "#1B6BC4",
                    E = "#A32B12", F = "#5C4FBF")

arm_cols_dark  <- c(A = "#2E9C81", C = "#D0762F", D = "#A277D6",
                    B2 = "#97941F", B3 = "#D44E92", B1 = "#4B92DB",
                    E = "#C9553C", F = "#8E82E8")

arm_cols <- arm_cols_light   # vignettes render on a light ground
