# Growth and Development - Solow, RCK, and OLG Models (South Korea)

Project applying the Solow, Ramsey-Cass-Koopmans (RCK), and overlapping generations (OLG) growth models to South Korea, using Penn World Table 11.0 data. Also runs some Barro-style convergence regressions on the same country over time.

## Data

Uses Penn World Table 11.0 (`pwt110.xlsx`), filtered to Republic of Korea, from roughly 1960 onward. You'll need to download PWT 11.0 yourself and point `read_excel()` at your own local path - the script currently has my own file path hardcoded near the top, so change that first.

Variables pulled from PWT:
- `rgdpe` - real GDP (expenditure side, PPP)
- `cn` - capital stock
- `emp` - employment (labour force proxy)
- `hc` - human capital index
- `csh_i` - investment share of GDP (savings rate)
- `irr` - real interest rate
- `delta` - depreciation rate
- `csh_c` - consumption share of GDP
- `rtfpna` - TFP growth
- `rconna` - real consumption
- `labsh` - labour share of output

## What's in the script

1. Load and filter the PWT data to Korea
2. Build per-worker terms (y, k, c) and growth rates (n, g)
3. Calibrate parameters: alpha (capital share), theta (inverse IES, set to 2), and back out rho and A
4. Solow model - actual vs break-even investment, steady state k* and y*
5. RCK model - Euler equation, k-dot and c-dot loci, phase diagrams
6. Barro convergence regressions (unconditional and conditional), run on Korea over time rather than a cross-country panel, so it's really testing mean reversion
7. OLG model - law of motion, transition function, steady state, and a dynamic inefficiency check (r vs n+g)

Plots for all of the above are built in as you go using ggplot2, styled with a WSJ-style theme from `ggthemes`.

## Packages needed

```r
install.packages(c("tidyverse", "readxl", "estimatr", "zoo", "writexl", "ggthemes"))
```

## Notes

This was for a growth and development course, mostly to get a feel for how these models actually look when you throw real data at them instead of just doing the algebra. The Barro regressions here are single-country over time, which isn't the textbook cross-country version, so don't read too much into the coefficients - it's more about seeing the mechanics work.

## Limitations

k\* (and y\*) are calculated once from the sample-average alpha, rho, and g, and then treated as a single fixed target for the whole time period. In reality g and alpha aren't constant over 60+ years, so the "true" steady state is probably moving around too, not sitting still. I think this was the main limitation with the approach here, though I'm not 100% sure that's the full picture - just flagging it as something to be aware of if you're using these plots to read too much into the convergence story.
