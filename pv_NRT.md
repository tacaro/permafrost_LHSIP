# Calculating partial pressures of gas

We know the conc of analyte in gas standard and we know the volume and the pressure.
With ideal gas law we can get to moles of gas

$$ 
PV = nRT
$$

GC measures moles, nothing else.

When drawing sample from samples: did not vent. We know the Volume and the moles.
We have moles div/volume to get concentration via ideal gas law. Use temp of room air.

$$
P = nRT/V
$$
```{r}
calculate_pressure <- function(
    n_moles,
    R = 8.401, # change me! add units!
    temp_C = 22, # change me!
    vol_mL) {
  pressure_atm = n_moles * R * temp_C / vol_mL
  return(pressure_atm)
}
```


Pressure want to use Boulder atmosphere (less than 1 atm).

For std curve: vol injected --> moles injected (n).
P = 0.8 atm
T = 20
R = 0.0821 L.atm.K-1.mol-1


$$
n = PV / RT
$$

```{r}
calculate_moles_gas <- function(
    pressure_atm = 1, # change me!
    volume_mL, 
    R = 8.401, # change me! add units!
    temp_C = 22 # change me!
    ) {
  n_moles = pressure_atm * volume_mL / R * temp_c
  return(n_moles)
}
```
