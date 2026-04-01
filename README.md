# econ_notes

Lecture Notes: Nested Logit Demand in Differentiated Products Markets

These notes synthesize the two lecture PDFs — Nested (1).pdf and Nested_2.pdf — covering logit demand, the IIA problem, nested logit, estimation, marginal cost recovery, merger simulation, and applications.



1. Starting Point: The Standard Logit Model

Setup

We have a market with $J$ differentiated products and one outside good (product 0, e.g. "don't buy"). Consumer $i$ gets utility from product $j$:

$$U_{ij} = \delta_j + \varepsilon_{ij}$$

where:





$\delta_j = X_j \beta + \alpha p_j + \xi_j$ is mean utility — the part common to all consumers





$X_j$: observed product chara/cteristics (range, speed, volume, etc.)



$p_j$: price



$\alpha < 0$: price coefficient (higher price → lower utility)



$\xi_j$: unobserved product quality (the structural error term)



$\varepsilon_{ij}$: idiosyncratic taste shock, drawn i.i.d. from a Type I Extreme Value (Gumbel) distribution

The outside good has utility $U_{i0} = \varepsilon_{i0}$ (normalized mean utility = 0).

Market Shares

Under the Gumbel distribution, the probability that consumer $i$ chooses product $j$ is:

$$s_j = \frac{\exp(\delta_j)}{1 + \sum_{k=1}^J \exp(\delta_k)}$$

The "1" in the denominator comes from the outside good ($\exp(0) = 1$).

Berry (1994) Inversion

A key insight: we can invert the share equation to get a linear estimating equation. Take logs:

$$\ln s_j = \delta_j - \ln\left(1 + \sum_k \exp(\delta_k)\right)$$

$$\ln s_0 = -\ln\left(1 + \sum_k \exp(\delta_k)\right)$$

Subtract:

$$\ln s_j - \ln s_0 = \delta_j = X_j \beta + \alpha p_j + \xi_j$$

This is a linear regression where the dependent variable $\ln s_j - \ln s_0$ is constructed from the data, and $\xi_j$ is the error term.

Endogeneity Problem

We cannot run OLS on this equation because price $p_j$ is correlated with $\xi_j$. Why? Firms observe $\xi_j$ (their product's unobserved quality) and set prices accordingly — high-quality products command higher prices. OLS would give a biased (less negative) estimate of $\alpha$.

We need instruments for price: variables that shift the supply side (and thus prices) but don't directly affect demand. Examples include cost shifters, characteristics of competing products, etc.

Marginal Cost Recovery (Single-Product Firms)

Under Bertrand-Nash competition, each single-product firm $j$ maximizes:

$$\pi_j = (p_j - c_j) \cdot s_j(p) \cdot M$$

where $M$ is market size and $c_j$ is marginal cost.

The first-order condition is:

$$s_j + (p_j - c_j) \frac{\partial s_j}{\partial p_j} = 0$$

For the logit:

$$\frac{\partial s_j}{\partial p_j} = \alpha \cdot s_j(1 - s_j)$$

Substituting:

$$s_j + (p_j - c_j) \cdot \alpha \cdot s_j(1 - s_j) = 0$$

$$1 + (p_j - c_j) \cdot \alpha (1 - s_j) = 0$$

$$\boxed{c_j = p_j + \frac{1}{\alpha(1 - s_j)}}$$

Since $\alpha < 0$, the second term is negative, so $c_j < p_j$ (positive markups), as expected.

The IIA Problem

The logit model has a well-known limitation: Independence of Irrelevant Alternatives (IIA). The ratio of any two products' shares depends only on their own characteristics:

$$\frac{s_j}{s_k} = \frac{\exp(\delta_j)}{\exp(\delta_k)} = \exp(\delta_j - \delta_k)$$

This ratio does not depend on any other product in the market. This implies:





Unrealistic substitution patterns: If a product is removed, its market share is redistributed proportionally to all other products (including the outside good). A removed sports car would lose equal proportional share to a minivan and to another sports car.



Cross-price elasticities depend only on the rival's market share, not on how similar two products are.

This is unrealistic. In practice, consumers substitute more toward products that are similar (e.g., same segment). The nested logit fixes this.



2. The Nested Logit Model

Motivation

The nested logit introduces correlation among products within the same group (nest). Products in the same segment (e.g., small jets, medium jets, large jets) are closer substitutes to each other than to products in different segments. This breaks IIA across groups while preserving it within groups.

Model Setup

Products are partitioned into $G$ nests (groups). In our application, nests are jet segments: $g \in S, M, L$. The outside good is in its own singleton nest.

Consumer $i$'s utility for product $j$ in nest $g$:

$$U_{ij} = \delta_j + \zeta_{ig} + (1 - \lambda)\varepsilon_{ij}$$

where:





$\delta_j = X_j \beta + \alpha p_j + \xi_j$: mean utility (same as before)



$\zeta_{ig}$: a nest-specific random shock — captures unobserved preference for segment $g$ (e.g., consumer $i$ prefers large jets)



$\varepsilon_{ij}$: idiosyncratic product-specific shock



$\lambda \in (0, 1]$: nesting parameter controlling within-nest correlation

Interpreting $\lambda$

The nesting parameter $\lambda$ is crucial:







Value of $\lambda$



Meaning





$\lambda = 1$



No correlation within nests → collapses to standard logit





$\lambda \to 0$



Perfect correlation within nests → consumers first choose a segment, then pick uniformly within it





$0 < \lambda < 1$



Products in the same nest are closer substitutes than products in different nests

Technically, the error term $\zeta_{ig} + (1-\lambda)\varepsilon_{ij}$ has a Generalized Extreme Value distribution that generates the nested logit probabilities. The key mathematical property: the joint distribution allows for correlation within nests while maintaining closed-form share expressions.

Share Formulas

The nested logit share has a multiplicative structure:

$$s_j = s_{j|g} \cdot s_g$$

Within-nest share (conditional on choosing nest $g$):

$$s_{j|g} = \frac{\exp(\delta_j / \lambda)}{D_g}$$

where $D_g = \sum_{k \in g} \exp(\delta_k / \lambda)$ is the inclusive value (or log-sum) of nest $g$.

Nest share (probability of choosing nest $g$):

$$s_g = \frac{D_g^\lambda}{1 + \sum_h D_h^\lambda}$$

The "1" again comes from the outside good.

Outside good share:

$$s_0 = \frac{1}{1 + \sum_h D_h^\lambda}$$

Intuition for $D_g$

The inclusive value $D_g = \sum_{k \in g} \exp(\delta_k/\lambda)$ summarizes how attractive nest $g$ is. More products, higher quality products, or lower prices all increase $D_g$, which increases $s_g$ — making the nest more popular.

The exponent $\lambda$ on $D_g$ in the nest share formula controls how sensitive nest choice is to the inclusive value. When $\lambda$ is small, the nest share is very sensitive to the quality of the best option within the nest.

Berry Inversion for Nested Logit

Following the same log-subtraction trick as before:

$$\ln s_j = \frac{\delta_j}{\lambda} + (\lambda - 1)\ln D_g - \ln(1 + \sum_h D_h^\lambda)$$

$$\ln s_0 = -\ln(1 + \sum_h D_h^\lambda)$$

Subtract:

$$\ln s_j - \ln s_0 = \frac{\delta_j}{\lambda} + (\lambda - 1)\ln D_g$$

Now notice that $\ln s_{j|g} = \delta_j/\lambda - \ln D_g$, so:

$$\ln s_j - \ln s_0 = \underbrace{\ln s_{j|g} + \ln D_g}{\delta_j/\lambda} + (\lambda - 1)\ln D_g = \ln s{j|g} + \lambda \ln D_g$$

We can also write this as:

$$\ln s_j - \ln s_0 = \delta_j + (1 - \lambda) \ln s_{j|g}$$

This is the key estimating equation:

$$\boxed{\ln s_j - \ln s_0 = X_j \beta + \alpha p_j + (1-\lambda) \ln s_{j|g} + \xi_j}$$

Defining $\sigma = 1 - \lambda$, this is:

$$\ln s_j - \ln s_0 = X_j \beta + \alpha p_j + \sigma \ln s_{j|g} + \xi_j$$

where:





LHS: log share ratio (computed from data)



$X_j \beta$: product characteristics



$\alpha p_j$: price effect (endogenous)



$\sigma \ln s_{j|g}$: within-nest share effect (endogenous)



$\xi_j$: unobserved quality (error term)

Two Endogenous Variables

Both price and $\ln s_{j|g}$ are endogenous:





Price: correlated with $\xi_j$ because firms set higher prices for higher unobserved quality



$\ln s_{j|g}$ (within-nest share): correlated with $\xi_j$ because a product with higher unobserved quality captures a larger share within its nest

We therefore need at least two instruments — one for each endogenous variable.



3. Estimation Strategy: 2SLS with BLP Instruments

The Instrument Problem

We need variables that:





Are correlated with price and within-nest share (relevance)



Are uncorrelated with unobserved product quality $\xi_j$ (exogeneity)

BLP Instruments (Berry, Levinsohn & Pakes, 1995)

The idea: the characteristics of competing products affect equilibrium prices and shares through competition but are plausibly exogenous to any single product's unobserved quality.

For product $j$ produced by firm $f(j)$ in market $t$:

$$Z_{jt} = \sum_{\substack{k \neq j  f(k) \neq f(j)}} X_{kt}$$

This is the sum of rival firms' product characteristics in the same market. We construct one instrument for each characteristic (range, volume, speed, power, etc.).

Why do these work?





Relevance: More (or better) rival products increase competitive pressure, which lowers equilibrium prices and shifts within-nest shares. So these instruments are correlated with the endogenous variables.



Exogeneity: Rival product characteristics are determined by rivals' technology and design choices, which are plausibly uncorrelated with product $j$'s unobserved quality $\xi_j$.

2SLS Procedure

Step 1 — First stage: Regress each endogenous variable (price, $\ln s_{j|g}$) on all exogenous regressors (constant, $X_j$, segment dummies) and the excluded BLP instruments. Get fitted values $\hat{p}j$ and $\widehat{\ln s{j|g}}$.

Step 2 — Second stage: Replace the endogenous variables with their first-stage fitted values and run OLS:

$$\ln s_j - \ln s_0 = X_j \beta + \alpha \hat{p}j + \sigma \widehat{\ln s{j|g}} + \text{error}$$

Equivalently, the 2SLS estimator in matrix form:

$$\hat{\beta}_{2SLS} = (X'P_Z X)^{-1} X'P_Z y$$

where $P_Z = Z(Z'Z)^{-1}Z'$ is the projection onto the instrument space.

Standard errors use the actual residuals $\hat{u} = y - X\hat{\beta}$ (not the second-stage residuals):

$$\hat{V}(\hat{\beta}) = \hat{\sigma}^2 (X' P_Z X)^{-1}, \quad \hat{\sigma}^2 = \frac{\hat{u}'\hat{u}}{n - K}$$

First-Stage Diagnostics

Always check:





Partial F-statistic for excluded instruments in each first-stage regression. Rule of thumb: $F > 10$ suggests instruments are not weak.



R-squared of first-stage regressions.

Weak instruments lead to biased 2SLS estimates (biased toward OLS).



4. Price Derivatives in the Nested Logit

The derivatives of shares with respect to prices are essential for computing marginal costs via the FOCs and for merger simulation.

Since $\delta_j = X_j\beta + \alpha p_j + \xi_j$, we have $\partial \delta_j / \partial p_j = \alpha$.

Own-Price Derivative

$$\frac{\partial s_j}{\partial p_j} = \alpha \cdot s_j \left[\frac{1}{\lambda}(1 - s_{j|g}) + s_{j|g}(1 - s_g)\right]$$

Since $\alpha < 0$ and the bracketed term is positive, own-price derivatives are negative (raising your price lowers your share).

Interpretation of the two terms inside the bracket:





$\frac{1}{\lambda}(1 - s_{j|g})$: within-nest substitution effect. When $\lambda$ is small, this term is large — a price increase causes big losses to other products in the same nest.



$s_{j|g}(1 - s_g)$: across-nest substitution. Some consumers leave the nest entirely.

Cross-Price Derivative: Same Nest

For products $j$ and $k$ in the same nest $g$ ($j \neq k$):

$$\frac{\partial s_k}{\partial p_j} = -\alpha \cdot s_k \cdot s_{j|g} \left(\frac{1}{\lambda} - 1 + s_g\right)$$

Since $\alpha < 0$, this is positive: raising $j$'s price increases $k$'s share (consumers substitute toward $k$).

The key term is $(1/\lambda - 1)$. When $\lambda < 1$, this is positive and potentially large, meaning within-nest cross-price effects are stronger than across-nest effects.

Cross-Price Derivative: Different Nests

For $j$ in nest $g$ and $k$ in nest $h \neq g$:

$$\frac{\partial s_k}{\partial p_j} = -\alpha \cdot s_k \cdot s_j$$

This is the same as in the standard logit — it depends only on the products' market shares, not on any within-nest terms.

Comparison: Within-Nest vs. Across-Nest Substitution

The ratio of within-nest to across-nest cross-price effects is approximately:

$$\frac{\text{within-nest effect}}{\text{across-nest effect}} \approx \frac{s_{j|g}(1/\lambda - 1 + s_g)}{s_j} = \frac{1/\lambda - 1 + s_g}{s_g}$$

When $\lambda$ is small (strong nesting), this ratio is large — products in the same segment are much closer substitutes. When $\lambda = 1$ (no nesting), the ratio equals 1, and we're back to the logit.

Verification: Collapse to Logit When $\lambda = 1$

Setting $\lambda = 1$:





Own-price: $\alpha \cdot s_j[(1)(1 - s_{j|g}) + s_{j|g}(1 - s_g)] = \alpha \cdot s_j(1 - s_j)$ ✓



Within-nest cross: $-\alpha \cdot s_k \cdot s_{j|g} \cdot s_g = -\alpha \cdot s_k \cdot s_j$ (equals the across-nest formula) ✓



5. Marginal Cost Recovery for Multi-Product Firms

Why Multi-Product Firms Matter

In the jets market, firms like Bombardier, Dassault, and Gulfstream sell products in multiple segments. When Bombardier sets the price for its medium jet, it accounts for how that price affects sales of its small and large jets too. This internalization of cannibalization leads to higher prices than if each product were priced independently.

The Multi-Product FOC System

Firm $f$ owns a set of products $\mathcal{F}_f$. It maximizes total profit:

$$\Pi_f = \sum_{j \in \mathcal{F}_f} (p_j - c_j) \cdot s_j(\mathbf{p}) \cdot M$$

The FOC for product $j \in \mathcal{F}_f$:

$$s_j + \sum_{k \in \mathcal{F}_f} (p_k - c_k) \frac{\partial s_k}{\partial p_j} = 0$$

The second term sums over all of firm $f$'s products — not just product $j$. If firm $f$ also owns product $k$, then the effect of $j$'s price on $k$'s sales is internalized.

Matrix Form

Define the ownership-weighted derivative matrix $\Omega$ (dimension $J \times J$):

$$\Omega_{jk} = \begin{cases} \frac{\partial s_k}{\partial p_j} & \text{if } j \text{ and } k \text{ are owned by the same firm}  0 & \text{otherwise} \end{cases}$$

The system of FOCs across all products is:

$$\mathbf{s}(\mathbf{p}) + \Omega(\mathbf{p}) \cdot (\mathbf{p} - \mathbf{c}) = \mathbf{0}$$

Solving for marginal costs:

$$\boxed{\mathbf{c} = \mathbf{p} + \Omega^{-1} \mathbf{s}}$$

Since $\Omega$ has negative diagonal entries (own-price effects) and positive off-diagonal entries (cross-price effects for same-firm products), $\Omega^{-1}\mathbf{s}$ is a negative vector, giving $\mathbf{c} < \mathbf{p}$ (positive markups).

Single-Product Firm Special Case

If firm $f$ owns only product $j$, the matrix $\Omega$ restricted to $j$ is just the scalar $\partial s_j / \partial p_j$. The FOC reduces to:

$$c_j = p_j + \frac{s_j}{\partial s_j / \partial p_j} = p_j + \frac{1}{\alpha[1/\lambda(1-s_{j|g}) + s_{j|g}(1-s_g)]}$$

For $\lambda = 1$ (logit), this further simplifies to $c_j = p_j + 1/[\alpha(1 - s_j)]$.

Diversion Ratios

The diversion ratio from product $j$ to product $k$ measures what fraction of $j$'s lost sales go to $k$ when $j$'s price increases:

$$D_{j \to k} = \frac{-\partial s_k / \partial p_j}{\partial s_j / \partial p_j}$$

In the nested logit:





Within-nest diversion is higher than across-nest diversion (products in the same segment are closer substitutes)



The logit's IIA implies equal diversion to all alternatives (proportional to share), which is unrealistic



The nested logit improves by creating asymmetric diversion: more diversion toward same-segment products

Diversion ratios are important in antitrust because they measure how directly two products compete. High diversion between merging firms' products implies larger price increases from a merger.



6. Merger Simulation

Framework

A merger changes the ownership structure. Products that were priced independently are now priced jointly. The merged firm internalizes cross-product effects, typically raising prices (absent cost savings).

Steps





Estimate demand (get $\alpha$, $\beta$, $\lambda$)



Recover marginal costs from pre-merger prices using the pre-merger ownership $\Omega^{pre}$



Construct post-merger ownership $\Omega^{post}$ (merging firms' products now internalized)



Solve for new equilibrium prices under $\Omega^{post}$ with same marginal costs



Compare pre- and post-merger prices, shares, and welfare

Solving for Post-Merger Equilibrium: Iterated Best Response

The equilibrium condition is:

$$\mathbf{p} = \mathbf{c} - \Omega(\mathbf{p})^{-1} \mathbf{s}(\mathbf{p})$$

This is a fixed-point equation. The iterated best response algorithm:





Start with pre-merger prices: $\mathbf{p}^{(0)} = \mathbf{p}^{pre}$



Compute shares and derivatives at $\mathbf{p}^{(n)}$



Update: $\mathbf{p}^{(n+1)} = \mathbf{c} - \Omega(\mathbf{p}^{(n)})^{-1} \mathbf{s}(\mathbf{p}^{(n)})$



Repeat until $\mathbf{p}^{(n+1)} - \mathbf{p}^{(n)} < \varepsilon$

This typically converges in a few dozen iterations.

Merger with Efficiency Gains

Mergers may also generate cost savings (economies of scale, elimination of duplicate operations). If marginal costs fall enough, the merger could lower prices despite increased market power.

To find the break-even efficiency: search for the cost reduction $\tau$ such that all post-merger prices equal or fall below pre-merger prices when applying $c_j^{post} = (1-\tau) c_j^{pre}$ for the merging firms' products.



7. Counterfactual: New Product Entry

Setup

To evaluate whether a firm should introduce a new product:





Define the new product: assign it average characteristics of the target segment, set $\xi = 0$ (no unobserved quality advantage), and assign marginal cost equal to the segment average



Expand the market: add the new product to the choice set



Solve for the new equilibrium: all firms (including the entrant) re-optimize prices



Compare profits: is the firm better off with the new product?

Key Economic Forces





Market expansion: the new product attracts some consumers from the outside good, growing total sales



Business stealing: the new product takes share from existing products, especially within the same nest



Cannibalization: if the entering firm already has a product in the target segment, the new product may steal from its own existing products



Strategic interaction: other firms respond by adjusting their prices (typically lowering them in the entry segment)

The optimal segment to enter balances these forces. Entering a segment where the firm has no existing products avoids cannibalization. Entering a segment with few competitors or high margins may yield the largest profit gain.



8. Two-Layer Nested Logit (Extension)

The lectures also cover a two-layer nested logit for richer substitution patterns. Instead of one level of nesting, products are grouped into sub-nests within nests.

Example (Verboven 1996, European car market):





Top nest: country of origin (domestic vs. foreign)



Sub-nest: market segment (subcompact, compact, midsize, luxury)

The estimation equation adds a second within-share term:

$$\ln s_j - \ln s_0 = X_j\beta + \alpha p_j + (1-\sigma_1)\ln s_{j|C} + (1-\sigma_2)\ln s_{C|G} + \xi_j$$

where $s_{j|C}$ is the within-sub-nest share, and $s_{C|G}$ is the sub-nest share within the top nest.

This requires three instruments (price, $\ln s_{j|C}$, $\ln s_{C|G}$ are all endogenous), so more BLP instruments are needed.



9. Applications from the Lectures

Einav (2007): Movie Seasonality

Uses nested logit to study the U.S. movie market. Movies are nested by genre. Key findings:





Within-genre substitution is much stronger than across-genre



Explains why Hollywood avoids releasing similar movies in the same week



Counterfactual: redistribution of release dates could increase industry revenue

Verboven (1996): European Car Market

Uses a two-layer nested logit to study car demand:





Products nested by domestic/foreign, then by segment



Finds domestic bias (consumers prefer home-country cars)



Uses the model to analyze trade policy counterfactuals



10. Summary of Key Formulas







Concept



Formula





Mean utility



$\delta_j = X_j\beta + \alpha p_j + \xi_j$





Within-nest share



$s_{j





Inclusive value



$D_g = \sum_{k \in g} \exp(\delta_k/\lambda)$





Nest share



$s_g = D_g^\lambda / (1 + \sum_h D_h^\lambda)$





Overall share



$s_j = s_{j





Berry inversion



$\ln s_j - \ln s_0 = X_j\beta + \alpha p_j + (1-\lambda)\ln s_{j





BLP instrument



$Z_{jt} = \sum_{k: f(k)\neq f(j)} X_{kt}$





Own-price deriv.



$\partial s_j/\partial p_j = \alpha s_j[1/\lambda(1-s_{j





Cross-price (same nest)



$\partial s_k/\partial p_j = -\alpha s_k s_{j





Cross-price (diff. nest)



$\partial s_k/\partial p_j = -\alpha s_k s_j$





Marginal cost recovery



$\mathbf{c} = \mathbf{p} + \Omega^{-1}\mathbf{s}$





Equilibrium fixed point



$\mathbf{p} = \mathbf{c} - \Omega(\mathbf{p})^{-1}\mathbf{s}(\mathbf{p})$



11. Glossary





IIA (Independence of Irrelevant Alternatives): Property of the logit where the ratio of two products' shares doesn't depend on other products. Implies unrealistic substitution patterns.



Nesting parameter ($\lambda$): Controls correlation within nests. $\lambda = 1$ gives the logit; smaller $\lambda$ means stronger within-nest substitution.



Inclusive value ($D_g$): Summarizes the attractiveness of a nest. Higher $D_g$ means the nest has better or more products.



Berry inversion: Technique to transform nonlinear share equations into a linear estimating equation by differencing with the outside good share.



BLP instruments: Sums of rival firms' product characteristics, used to instrument for endogenous price and within-share.



Diversion ratio: Fraction of product $j$'s lost sales captured by product $k$ when $j$'s price rises.



Iterated best response: Algorithm to find Nash equilibrium by repeatedly having each firm optimize given others' prices.



Cannibalization: When a firm's new product steals sales from its own existing products.



Outside good: The option of not purchasing any product in the market. Its share pins down the overall level of demand.

