# CC 03 - Electricity and Magnetism

## Practicals [file](CC03_Electricity_and_Magnetism_Practicals.ipynb)
Practicals:
* RC Circuit. (Investigation of Capacitance)
* LR Circuit. (Investigation of Inductance)

## LR circuit [file](LR_Circuit_SKP.ipynb)

LR circuit with DC source:

* For growth of current,
$$ iR + L\frac{di}{dt} = V $$
$$ \implies \frac{di}{dt} = -\frac{iR - V}{L} $$

* For decay of current,
$$ iR + L\frac{di}{dt} = 0 $$
$$\implies \frac{di}{dt} = -\frac{iR}{L} $$

## RC circuit [file](RC_Circuit_SKP.ipynb)

RC circuit with DC source:

* For charging,
$$ iR + \frac{q}{C} = V $$
$$\implies \frac{dq}{dt} = -\frac{q - CV}{RC} $$

* For discharging,
$$ iR + \frac{q}{C} = 0 $$
$$ \implies \frac{dq}{dt} = -\frac{q}{RC} $$

RC circuit with AC source:

* For charging,
$$ iR + \frac{q}{C} = V_0\sin(\omega t) $$
$$ \implies \frac{dq}{dt} = -\frac{q - CV_0\sin(\omega t)}{RC} $$

## LCR circuit [file](LCR_circuit_SKP.ipynb)

LCR circuit in DC supply:

* Charging of capacitor,
$$ L\frac{di}{dt} + \frac{q}{C} + iR = V $$
$$ \implies q''+2bq'+\omega_0^2q = \frac{V}{L}  $$
where $2b=R/L, \omega_0^2=1/LC$.

* Discharging of capacitor, 
$$ L\frac{di}{dt} + \frac{q}{C} + iR = 0 $$
$$ \implies q''+2bq'+\omega_0^2q = 0  $$
where $2b=R/L, \omega_0^2=1/LC$.

## LC circuit [file](LC_Circuit_SKP.ipynb)

LC circuit with DC source:

$$ L\frac{di}{dt} + \frac{q}{C} = 0 $$
$$ \implies \frac{d^2q}{dt^2} + \frac{q}{LC} = 0 $$
The current is
$$ i = \frac{dq}{dt} $$

