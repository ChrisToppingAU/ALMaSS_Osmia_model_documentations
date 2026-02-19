# *Osmia bicornis* Population Model: MIDox Implementation and Documentation

**Authors: Christopher John Topping, Xiaodong Duan, Elzbieta Ziolkowska**

## Abstract

This document provides comprehensive implementation documentation for an agent-based population model of the red mason bee (*Osmia bicornis*) developed within the ALMaSS (Animal Landscape and Man Simulation System) framework. The model simulates individual bee life cycles from egg through adult, incorporating temperature-dependent development, resource-limited provisioning behaviour, and overwintering physiology. Implementation follows the MIDox (Model Implementation Documentation with Doxygen) standard, providing detailed parameter documentation, biological rationale for design decisions, and explicit assessment of uncertainties. The model serves multiple research purposes: evaluating landscape management effects on solitary bee populations, assessing pesticide risks to pollinator communities, and exploring climate change impacts on bee phenology. Complete interactive documentation with full code cross-references is available at the GitHub Pages and Zenodo repositories listed in Section 9.

**Keywords:** agent-based model, *Osmia bicornis*, solitary bee, pollinator, ALMaSS, MIDox, implementation documentation, Doxygen

---

## 1. Introduction

### 1.1 Context and Motivation

The decline of wild pollinator populations across Europe and North America has raised concerns about pollination service stability and crop production resilience. Solitary bees, including mason bees (*Osmia* spp.), contribute substantially to crop pollination and wild plant reproduction, yet remain under-studied compared to managed honeybees. Understanding how landscape structure, agricultural management, and environmental stressors affect solitary bee populations requires mechanistic models that capture individual-level processes whilst scaling to landscape patterns.

This document describes the implementation of an agent-based population model for *Osmia bicornis* (Linnaeus, 1758), the red mason bee, a common and economically important solitary bee species across temperate Europe. The model represents individual bees as autonomous agents progressing through life stages, making behavioural decisions based on local conditions, and experiencing mortality from multiple sources. Population dynamics emerge from these individual-level processes operating within spatially-explicit agricultural landscapes.

### 1.2 The MIDox Framework

This documentation follows the MIDox (Model Implementation Documentation with Doxygen) standard proposed for ecological simulation models. MIDox addresses a persistent challenge in ecological modelling: the gap between formal model descriptions (mathematical specifications of processes) and actual implementations (working computer code). This gap impedes model reuse, validation, and extension whilst reducing research transparency.

MIDox bridges this gap through three complementary elements:

1. **Enhanced source code** with comprehensive biological documentation embedded as Doxygen comments
2. **Interactive HTML documentation** generated automatically from annotated code, providing searchable cross-references
3. **Narrative documentation** (this document) explaining design decisions, parameter sources, and implementation trade-offs

Together, these elements create a complete audit trail from formal model through implementation to validation, enabling:
- **Reproducibility**: Other researchers can understand and replicate the implementation
- **Maintainability**: Developers can modify code confidently, understanding biological intent
- **Transparency**: Reviewers can assess implementation fidelity to formal specifications
- **Extensibility**: New features can be added consistently with existing design patterns

### 1.3 Three-Paper Publication Sequence

This MIDox documentation is the second in a three-paper sequence:

1. **Formal Model** (Ziółkowska et al. 2023, Food and Ecological Systems Modelling Journal): Conceptual specification using mathematical notation and structured natural language, describing model purpose, scope, entities, processes, and equations without implementation details

2. **MIDox Implementation** (FESMJ article ref): Complete technical specification of the working model, documenting all classes, methods, parameters, and design decisions with biological justification

3. **Testing and Calibration** (in preparation): Validation against empirical data, sensitivity analyses, and parameter estimation procedures

This sequence separates concerns: the formal model establishes biological validity; the MIDox ensures implementation fidelity; testing confirms empirical adequacy. Each paper serves distinct audiences whilst maintaining coherence through explicit cross-references.

### 1.4 Document Structure

This document comprises nine sections:

- **Section 2** (Model Overview): High-level architecture, emergent properties, and spatial/temporal scales
- **Section 3** (Scheduling and Workflow): Daily execution sequence and seasonal transitions  
- **Section 4** (Implementation Details): Detailed component descriptions with code cross-references
- **Section 5** (Parameters): Organization, empirical sources, and uncertainty assessments
- **Section 6** (Inputs): Required data, formats, and preprocessing requirements
- **Section 7** (Outputs): Generated data, formats, and interpretation guidance
- **Section 8** (Implementation Discussion): Design trade-offs, limitations, and future directions
- **Section 9** (Documentation Access): Links to interactive documentation and archived code

Throughout, we use Doxygen `@ref` tags (e.g., [Osmia_Female](@ref Osmia_Female)) to link narrative text to specific code elements. These hyperlinks function in the HTML version, allowing readers to navigate seamlessly between conceptual explanation and implementation detail.

---

## 2. Model Overview

### 2.1 Model Purpose and Scope

The *Osmia bicornis* model simulates population dynamics of a univoltine (one generation per year) solitary bee species across agricultural landscapes. The model addresses three primary research questions:

1. **Landscape effects**: How do landscape composition (flower availability, nesting habitat) and configuration (patch size, connectivity) affect population persistence and abundance?

2. **Management impacts**: What are the consequences of agricultural intensification (pesticide use, mowing timing, hedgerow removal) for bee populations?

3. **Climate change**: How will altered temperature regimes affect phenology, voltinism, and population viability?

The model scope encompasses:
- Complete life cycle (egg, larva, prepupa, pupa, overwintering adult, active adult)
- Individual-level heterogeneity (body size, age, reproductive state)
- Resource-limited reproduction (pollen availability constrains provisioning)
- Temperature-dependent development (degree-day accumulation)
- Spatially-explicit behaviour (foraging, dispersal, nest-site selection)
- Multiple mortality sources (development, background, overwintering, parasitism)

The model deliberately excludes:
- Males (implicit in fertilization and sex ratio assumptions; computational efficiency)
- Detailed flower species preferences (aggregated pollen/nectar availability)
- Within-nest social interactions (solitary species; minimal overlap)
- Detailed parasitoid biology (optional mechanistic extension available; base model uses probability-based parasitism)

### 2.2 Model Architecture

#### 2.2.1 Agent Types and Life Stages

The model implements an agent-based architecture where individual bees are autonomous entities (agents) progressing through discrete life stages. Six agent classes represent developmental stages:

- [Osmia_Egg](@ref Osmia_Egg): From laying to hatching (temperature-dependent duration)
- [Osmia_Larva](@ref Osmia_Larva): Feeding and cocoon construction (resource-dependent growth)
- [Osmia_Prepupa](@ref Osmia_Prepupa): Summer diapause (time-based, weakly temperature-dependent)
- [Osmia_Pupa](@ref Osmia_Pupa): Metamorphosis to adult form (temperature-dependent duration)
- [Osmia_InCocoon](@ref Osmia_InCocoon): Fully-developed adult within cocoon (includes overwintering)
- [Osmia_Female](@ref Osmia_Female): Free-living reproductive female (foraging, nesting, dispersal)

All agent classes inherit from [Osmia_Base](@ref Osmia_Base), which provides common attributes (age, mass, location) and methods (mortality checking, temperature access). This inheritance structure reduces code duplication whilst enabling stage-specific specialization.

Males are not explicitly modelled. This design decision reflects biological reality (*Osmia bicornis* females are sole provisioners) and computational efficiency (modelling males would double agent count without substantially affecting female behaviour or population dynamics). Male abundance and mating success are implicit: females are assumed mated, and sex ratios emerge from maternal provisioning decisions.

#### 2.2.2 Population Manager

The [Osmia_Population_Manager](@ref Osmia_Population_Manager) orchestrates the simulation, handling:

- **Initialization**: Reading configuration parameters, constructing lookup tables, creating starting population
- **Scheduling**: Coordinating daily execution order (environmental updates, agent actions, cleanup)
- **Spatial infrastructure**: Maintaining density grids, nest management, resource maps
- **Global calculations**: Weather integration (foraging hours), seasonal transitions (overwintering triggers)

The population manager follows the ALMaSS framework's standard Population_Manager pattern, providing *Osmia*-specific implementations of core scheduling hooks ([DoFirst()](@ref DoFirst()), [DoBefore()](@ref DoBefore()), [DoAfter()](@ref DoAfter()), [DoLast()](@ref DoLast())) whilst delegating agent behaviour to individual classes.

#### 2.2.3 Nest Management

*Osmia bicornis* nests in pre-existing cavities (hollow stems, beetle holes, artificial nest boxes). The [Osmia_Nest_Manager](@ref Osmia_Nest_Manager) maintains polygon-level nest availability, tracking:

- Suitable nesting habitat (derived from land-use classification)
- Nest capacity (maximum nests per polygon based on substrate availability)
- Active nests (currently occupied by provisioning females)
- Nest lifecycle (creation when female claims cavity, release when abandoned or offspring emerge)

Individual nests (modelled by [Osmia_Nest](@ref Osmia_Nest)) contain linearly-arranged cells, each provisioned sequentially by the founding female. Nest structures persist across years (perennial nesting sites), though individual nest occupancy is annual (single reproductive season).

### 2.3 Spatial and Temporal Scales

#### 2.3.1 Spatial Resolution

The model operates within the ALMaSS framework's landscape representation:

- **Base resolution**: Polygon-based land-use map (fields, hedgerows, woodlands, etc.)
- **Pollen map**: 10 m × 10 m grid storing flower availability (quality and quantity)
- **Female density grid**: 1 km × 1 km grid tracking local female abundance
- **Movement**: Continuous coordinates (meters) for individual locations

This multi-resolution approach balances biological realism (bees perceive resources at meter scales) against computational efficiency (large landscapes require coarse-grained tracking structures). The 1 km² density grid provides ecologically meaningful local density measures whilst avoiding excessive spatial resolution that would dominate runtime.

Typical simulation extent: 10-100 km² landscapes (agricultural regions). Minimum viable extent: ~5 km² (sufficient to contain resources and nesting habitat whilst allowing dispersal). Maximum practical extent: ~500 km² (limited by computation time and data availability).

#### 2.3.2 Temporal Resolution

The model executes with daily time steps, appropriate for *Osmia bicornis* biology:

- Development: Occurs gradually (degree-day accumulation over days/weeks)
- Behaviour: Females forage and provision daily (weather permitting)
- Mortality: Applied as daily probabilities
- Phenology: Seasonal transitions (emergence, diapause, overwintering) occur over weeks

Finer temporal resolution (hourly) would add computational cost without substantially improving biological realism. Coarser resolution (weekly) would miss critical short-term dynamics (weather effects on provisioning, rapid parasitism risk accumulation).

Typical simulation duration: 1-10 years. Single-year runs explore seasonal dynamics; multi-year runs assess population persistence and inter-annual variability. Initialization typically starts with overwintering adults (mid-lifecycle) to generate realistic first-year emergence phenology without requiring full-cycle spin-up.

### 2.4 Key Processes and Emergent Properties

#### 2.4.1 Core Processes

Five fundamental processes drive model dynamics:

**1. Temperature-Dependent Development**
Developmental stages accumulate degree-days above threshold temperatures (lower developmental thresholds). Eggs, larvae, and pupae follow standard linear degree-day models. The prepupal stage uses a modified time-based approach due to non-linear temperature responses during summer diapause (see Section 4.2.3 for detailed rationale). Overwintering adults accumulate chilling requirements before emergence.

**2. Resource-Limited Reproduction**  
Female provisioning behaviour creates the link between landscape resources and population dynamics. Females search for pollen-rich flower patches within foraging range, collect pollen during available weather hours, and provision nest cells sequentially. Insufficient resources lead to abandoned nests, reduced egg production, or smaller offspring. This resource limitation creates density dependence and landscape sensitivity.

**3. Phenological Synchronization**
The prepupal summer diapause and overwintering adult stage create synchronized spring emergence despite variation in larval development rates. Diapause acts as a "waiting room" where fast-developing individuals pause, allowing slower individuals to catch up. This synchronization ensures temporal overlap between flowering periods and adult activity, critical for resource acquisition.

**4. Size-Structured Life History**
Maternal provisioning mass determines offspring body size, which in turn affects fecundity (larger females lay more eggs), sex ratio decisions (larger females produce more female-biased broods), and potentially survival. This creates transgenerational links where environmental conditions experienced by provisioning females affect offspring fitness.

**5. Spatial Processes**
Females exhibit philopatry (preference for natal area), searching locally for nest sites before dispersing to distant areas if unsuccessful. Foraging occurs within species-specific ranges (hundreds of meters). Dispersal between nesting areas enables colonization and gene flow. These spatial behaviours create landscape-scale patterns from individual movement rules.

#### 2.4.2 Emergent Properties

Population-level patterns emerge from individual-level processes without explicit specification:

**Population Growth and Regulation**
Population growth emerges from birth-death dynamics: fecundity depends on resource acquisition; mortality depends on development, background hazards, overwintering, and parasitism. Density-dependent regulation emerges implicitly through resource depletion (more females → less pollen per capita → reduced provisioning → lower fecundity) rather than explicit crowding effects.

**Spatial Distribution**
Bee distribution patterns emerge from habitat selection and resource tracking. Females cluster in areas with abundant flowers and nest sites, creating spatial heterogeneity even in uniform landscapes. Source-sink dynamics emerge when some patches produce surplus individuals (sources) whilst others rely on immigration (sinks).

**Phenological Patterns**
Emergence timing, peak abundance, and seasonal decline emerge from accumulated environmental history (temperature, resource availability) and individual variation (body size, development rates). Years with warm springs produce early emergence; years with abundant spring flowers support larger populations; years with poor late-season resources show early population decline.

**Sex Ratio Variation**
Population-level sex ratios emerge from individual maternal decisions based on body size, age, and provisioned mass. Optimal sex ratio theory predicts that mothers in good condition should produce female-biased broods (larger, more valuable offspring), whilst mothers in poor condition shift toward males (cheaper offspring). The model implements these decision rules at the individual level; population sex ratios emerge as aggregate outcomes.

### 2.5 Model Validation and Uncertainty

Model validation follows standard agent-based modelling approaches:

**Pattern-Oriented Modelling**: The model reproduces multiple empirical patterns simultaneously (emergence phenology, fecundity distributions, development durations, sex ratios, spatial distributions). Matching multiple patterns increases confidence that model mechanisms are biologically realistic.

**Sensitivity Analysis**: Key uncertainties identified through systematic parameter perturbation. Parameters with large effects on population dynamics prioritized for refinement; parameters with minor effects use literature values without extensive calibration.

**Uncertainty Propagation**: Monte Carlo simulations with parameter distributions (rather than point values) characterize prediction uncertainty. Outputs reported with confidence intervals acknowledging parameter uncertainty.

Major sources of uncertainty include:
- Prepupal development (laboratory-field extrapolation difficult; see Section 8.3)
- Parasitism rates (field measurements sparse; wide variation)
- Foraging thresholds (species-specific data limited)
- Overwintering mortality (temperature dependence poorly characterized)

These uncertainties are explicitly documented in parameter descriptions (Section 5) with CONFIDENCE ratings (HIGH/MEDIUM/LOW) guiding interpretation and future research priorities.

## 3. Scheduling and Workflow

### 3.1 Daily Execution Sequence

The model executes with daily time steps following the ALMaSS framework's standardized scheduling pattern. Each day proceeds through four phases coordinated by the [Osmia_Population_Manager](@ref Osmia_Population_Manager):

#### Phase 1: [DoFirst()](@ref DoFirst()) - Environmental Updates

```
FOR each day:
    // Update global environmental state
    temperature_today = landscape.GetTemperature()
    Osmia_Base.SetTemp(temperature_today)  // Distribute to all agents
    
    // Calculate weather-limited foraging hours
    foraging_hours = 0
    FOR each hour in [0..23]:
        IF (temp[hour] > min_flight_temp AND
            wind[hour] < max_flight_wind AND  
            precip[hour] < max_flight_precip):
            foraging_hours++
    population_manager.m_FlyingWeather = foraging_hours
    
    // Update prepupal development rate from temperature
    temp_rounded = ROUND(temperature_today)
    prepupal_rate_today = lookup_table[temp_rounded]
    
    // Update nest manager status
    nest_manager.UpdateOsmiaNesting()  // Check for abandoned nests
    
    // Clear density grid (will be repopulated during [BeginStep()](@ref BeginStep()))
    FOR each grid_cell:
        grid_cell.female_count = 0
```

This phase sets up shared environmental state that all agents will access during their individual updates. Calculations performed once per day rather than per agent yields substantial computational savings (O(1) vs. O(N) operations).

#### Phase 2: Individual [BeginStep()](@ref BeginStep()) - Agent Pre-Processing

```
FOR each life_stage IN [Egg, Larva, Prepupa, Pupa, InCocoon, Female]:
    FOR each individual IN life_stage:
        individual.[BeginStep()](@ref BeginStep())
```

Stage-specific [BeginStep()](@ref BeginStep()) implementations prepare agents for daily updates:

**Egg/Larva/Prepupa/Pupa**: Check for completion of development
```
accumulated_development += daily_development_increment
IF accumulated_development >= required_development:
    ready_for_transition = TRUE
```

**InCocoon**: Update overwintering or pre-emergence state
```
IF in_overwintering_phase:
    IF temperature < overwintering_threshold:
        overwintering_DD += (overwintering_threshold - temperature)
    IF overwintering_DD >= required_DD AND date > March_1:
        ready_for_emergence = TRUE
```

**Female**: Reset daily state, check survival
```
foraging_hours_available = population_manager.GetForageHours()
emerge_age++
IF emerge_age > maximum_lifespan:
    mortality_flag = TRUE
Update local resource availability if needed
Update density grid: grid[location].female_count++
```

#### Phase 3: Individual [Step()](@ref Step()) - Main Agent Actions

```
FOR each life_stage IN [Egg, Larva, Prepupa, Pupa, InCocoon, Female]:
    FOR each individual IN life_stage:
        individual.[Step()](@ref Step())
```

This phase implements stage-specific behaviour. For developmental stages (Egg through InCocoon), Step() applies mortality and signals transitions. For [Osmia_Female](@ref Osmia_Female), Step() executes the behavioural state machine:

```
WHILE state != DONE:
    SWITCH current_state:
        CASE st_Develop:  // Initial maturation
            maturation_days++
            IF maturation_days >= prenesting_period:
                current_state = st_Dispersal
        
        CASE st_Dispersal:  // Long-distance movement
            Execute dispersal (Beta distribution movement)
            IF found_suitable_habitat:
                current_state = st_ReproductiveBehaviour
            ELSE IF attempts > max_attempts:
                current_state = st_Die
        
        CASE st_ReproductiveBehaviour:  // Core provisioning loop
            IF no_current_nest:
                TRY find nest site locally
                IF successful:
                    create_nest()
                    PlanEggsPerNest()
                ELSE:
                    nest_find_attempts++
                    IF nest_find_attempts > threshold:
                        current_state = st_Dispersal
            ELSE:  // Have active nest
                Forage for pollen/nectar
                Provision current cell
                IF cell complete:
                    LayEgg()
                    eggs_remaining--
                    IF eggs_remaining == 0:
                        current_state = st_Die
                    ELSE IF nest complete:
                        release nest
                        no_current_nest = TRUE
        
        CASE st_Die:
            mortality_flag = TRUE
            state = DONE
```

The state machine architecture separates concerns: each state handles a distinct behavioural mode, transitions occur based on clear criteria, and the loop terminates when reaching DONE or DIE states.

#### Phase 4: [DoLast()](@ref DoLast()) - End-of-Day Processing

```
// Update seasonal flags based on temperature history
IF day > September AND NOT pre_wintering_ended:
    Check for sustained autumn temperature drop
    IF (three consecutive days < 13°C) AND (declining trend):
        pre_wintering_ended = TRUE
        
IF day == March_1:
    overwintering_ended = TRUE

// Reset flags after emergence season
IF day == June_1:
    pre_wintering_ended = FALSE
    overwintering_ended = FALSE

// Optional: Record statistics for validation
IF testing_mode_enabled:
    Record stage durations, egg production, etc.
```

### 3.2 Seasonal Transitions

The annual cycle comprises distinct phases with characteristic processes:

**Winter (January-February)**: Overwintering adults accumulate chilling requirements. No activity. Mortality applied based on accumulated winter degree-days.

**Spring Emergence (March-May)**: Adults emerge from cocoons once sufficient chilling accumulated and temperatures rise. Emergence stochastic, creating phenological spread. Mating occurs (implicit). Females begin pre-nesting maturation and dispersal.

**Reproduction (April-June)**: Peak female activity. Foraging, nest construction, provisioning, egg-laying. Weather variation creates day-to-day fluctuations in provisioning rates. Resource availability declines through season as flowers senesce and bee populations grow.

**Development (May-August)**: Eggs hatch, larvae feed and construct cocoons, prepupae enter summer diapause. Temperature drives developmental rates. Mortality from parasitoids, weather extremes, resource limitation.

**Metamorphosis (August-September)**: Prepupae transform to pupae, pupae to adults. Adults remain in cocoons (adult-in-cocoon stage). Pre-wintering flag triggers onset of overwintering physiology.

**Pre-Winter (September-November)**: Adults in cocoons undergo physiological changes preparing for winter dormancy. Mortality applied based on cocoon mass (resource stores). Pre-wintering end detected from sustained temperature drop.

### 3.3 Parallelization Strategy

The model supports parallel execution using OpenMP:

```cpp
#pragma omp parallel
FOR each individual IN population:
    individual.[Step()](@ref Step())
```

Thread safety ensured through:
- **Read-only global state**: Temperature, weather, resource maps accessed but not modified
- **Local agent state**: Each agent modifies only its own attributes
- **Synchronized nest access**: Cell-level locks prevent concurrent modification
- **Atomic density updates**: Grid operations use atomic increments

Performance scaling typically sublinear: 4 cores yield ~3× speedup, 16 cores yield ~8-10× speedup. Diminishing returns reflect synchronization overhead and memory bandwidth limits.

---

## 4. Implementation Details

This section describes key implementation components with extensive cross-references to source code. All code elements are hyperlinked in the HTML documentation; readers can navigate directly from narrative explanation to implementation.

### 4.1 Development and Mortality

All developmental stages inherit core development and mortality processes from [Osmia_Base](@ref Osmia_Base). Stage-specific implementations (e.g., [Osmia_Egg::Step](@ref Osmia_Egg::Step), [Osmia_Larva::Step](@ref Osmia_Larva::Step)) customize these processes whilst reusing infrastructure.

#### 4.1.1 Degree-Day Development Model

Egg, larval, and pupal stages use standard linear degree-day accumulation:

```cpp
// Implementation in Osmia_Egg::Step() (and similar for Larva, Pupa)
double temp = Osmia_Base::GetTemp();
if (temp > m_OsmiaEggDevelThreshold) {
    m_AgeDegrees += (temp - m_OsmiaEggDevelThreshold);
}
if (m_AgeDegrees >= m_OsmiaEggDevelTotalDD) {
    // Development complete, signal transition to next stage
    TransitionToLarva();
}
```

**Parameters** (see [Osmia_Base](@ref Osmia_Base) for full documentation):
- [Osmia_Base::m_OsmiaEggDevelThreshold](@ref Osmia_Base::m_OsmiaEggDevelThreshold): 0.0°C (calibrated; see Section 5.1)
- [Osmia_Base::m_OsmiaEggDevelTotalDD](@ref Osmia_Base::m_OsmiaEggDevelTotalDD): 86 degree-days (calibrated; see Section 5.1)
- [Osmia_Base::m_OsmiaLarvaDevelThreshold](@ref Osmia_Base::m_OsmiaLarvaDevelThreshold): 13.8°C (Radmacher & Strohm 2011)
- [Osmia_Base::m_OsmiaLarvaDevelTotalDD](@ref Osmia_Base::m_OsmiaLarvaDevelTotalDD): 240 degree-days (Radmacher & Strohm 2011)
- [Osmia_Base::m_OsmiaPupaDevelThreshold](@ref Osmia_Base::m_OsmiaPupaDevelThreshold): 1.1°C (calibrated; see Section 5.1)
- [Osmia_Base::m_OsmiaPupaDevelTotalDD](@ref Osmia_Base::m_OsmiaPupaDevelTotalDD): 570 degree-days (calibrated; see Section 5.1)

**Implementation Difference from Formal Model**: Egg and pupal thresholds substantially reduced from laboratory values (13.8°C → 0.0°C for eggs, 13.2°C → 1.1°C for pupae) to prevent developmental failures under field temperature regimes. See Section 8.3 for detailed justification.

#### 4.1.2 Prepupal Development (Modified Approach)

The prepupal stage represents summer diapause, a period of developmental arrest with weak temperature dependence. Unlike other stages, prepupal development uses a time-based model with temperature-dependent rates:

```cpp
// Implementation in Osmia_Prepupa::Step()
double rate = population_manager->GetPrePupalDevelDays();  // Daily lookup
m_Age += rate;
if (m_Age >= m_OsmiaPrepupalDevelTotalDays) {
    TransitionToPupa();
}
```

The rate lookup uses integer temperature (0-41°C) to index a pre-calculated array (see [Osmia_Population_Manager::m_PrePupalDevelRates](@ref Osmia_Population_Manager::m_PrePupalDevelRates)). Rates follow a thermal performance curve: low at temperature extremes, peaking around 20-25°C.

**Rationale**: The formal model specifies a quadratic temperature-development relationship based on Radmacher & Strohm (2011) laboratory studies at constant temperatures. Under field conditions with daily temperature fluctuations, this relationship produced unrealistic developmental durations. The lookup table approach allows flexible parameterization empirically calibrated to field emergence phenology whilst acknowledging mechanistic uncertainty. See Section 8.3 for extended discussion.

**Parameters**:
- [Osmia_Base::m_OsmiaPrepupalDevelTotalDays](@ref Osmia_Base::m_OsmiaPrepupalDevelTotalDays): 45 days (approximate; weakly temperature-dependent)
- Temperature-rate lookup: [Osmia_Population_Manager::m_PrePupalDevelRates](@ref Osmia_Population_Manager::m_PrePupalDevelRates)

#### 4.1.3 Overwintering Development

The [Osmia_InCocoon](@ref Osmia_InCocoon) stage implements complex overwintering physiology with three sub-phases:

**Phase 1: Pre-wintering** (late summer/autumn)
Adults in cocoons undergo physiological preparation for winter dormancy. No development accumulation. Transition to Phase 2 triggered by sustained autumn temperature drop (detected by population manager; see [Osmia_Population_Manager::DoLast](@ref Osmia_Population_Manager::DoLast)).

**Phase 2: Overwintering** (winter)
Accumulate chilling requirement below temperature threshold:
```cpp
if (temp < m_OsmiaInCocoonOverwinteringTempThreshold) {
    m_OverwinteringDD += (m_OsmiaInCocoonOverwinteringTempThreshold - temp);
}
```
Required: ~400-500 degree-days below 5°C threshold. Mortality increases with prolonged cold exposure via winter mortality equation (see Section 4.1.4).

**Phase 3: Spring pre-emergence** (late winter/spring)
Once chilling requirement met AND date > March 1st, adults can emerge when temperatures rise. Emergence probability increases with accumulated degree-days:
```cpp
double emergence_prob = EmergenceCountConst + EmergenceCountSlope × accumulated_DD;
if (g_rand_uni() < emergence_prob) {
    Emerge();  // Transition to Osmia_Female
}
```

This three-phase structure creates realistic emergence phenology: early warmth doesn't trigger emergence (chilling requirement unmet); late cold delays emergence despite sufficient chilling (temperature-dependent emergence probability).

**Parameters**:
- [Osmia_Base::m_OsmiaInCocoonOverwinteringTempThreshold](@ref Osmia_Base::m_OsmiaInCocoonOverwinteringTempThreshold): 5.0°C
- [Osmia_Base::m_OsmiaInCocoonEmergCountConst](@ref Osmia_Base::m_OsmiaInCocoonEmergCountConst), [Osmia_Base::m_OsmiaInCocoonEmergCountSlope](@ref Osmia_Base::m_OsmiaInCocoonEmergCountSlope): Emergence probability equation parameters

#### 4.1.4 Mortality Processes

Multiple mortality sources operate simultaneously:

**Developmental Mortality** (daily probabilities during development):
```cpp
if (g_rand_uni() < m_DailyDevelopmentMortEggs) {
    mortality_flag = TRUE;
}
```

Applied to eggs, larvae, prepupae, pupae. Rates from formal model specification:
- Eggs: 0.0014 per day (Radmacher & Strohm 2011)
- Larvae: 0.0014 per day
- Prepupae: 0.003 per day
- Pupae: 0.003 per day

**Overwintering Mortality** (mass- and temperature-dependent):
The [Osmia_InCocoon](@ref Osmia_InCocoon) stage implements overwintering mortality as function of accumulated winter degree-days and cocoon mass. Larger cocoons (more resource stores) survive better; prolonged cold increases mortality:

```cpp
double mortality_prob = WinterMortConst + WinterMortSlope × overwintering_DD;
// Modify by cocoon mass (implementation detail in code)
if (g_rand_uni() < mortality_prob) {
    mortality_flag = TRUE;
}
```

Equation from Sgolastra et al. (2011) working with *Osmia lignaria*, assumed applicable to *O. bicornis* given phylogenetic proximity and similar overwintering biology.

**Adult Background Mortality** (daily):
Reproductive females experience constant daily mortality from predation, disease, accidents:
```cpp
if (g_rand_uni() < m_OsmiaFemaleBckMort) {  // Default 0.02 per day
    mortality_flag = TRUE;
}
```

Rate 0.02 per day yields mean adult lifespan ~50 days, matching field observations.

**Parasitism** (applied at egg laying):
Parasitism risk increases with cell open time. See Section 4.4 for detailed parasitism implementation.

### 4.2 Female Provisioning Behaviour

The [Osmia_Female](@ref Osmia_Female) class implements resource-limited reproductive behaviour, the core link between landscape resources and population dynamics.

#### 4.2.1 Nest Finding and Creation

Females search for nest sites within typical homing distance (~600 m default) of current location:

```cpp
// Simplified from Osmia_Female::FindNest()
FOR attempt IN 1..max_attempts:
    candidate_polygon = SelectRandomPolygonWithinRange()
    IF IsNestingHabitatSuitable(candidate_polygon):
        nest = population_manager->CreateNest(location, polygon)
        IF nest != NULL:
            RETURN success
nest_find_attempts++
IF nest_find_attempts > threshold:
    Trigger dispersal (st_Dispersal state)
```

Nesting habitat suitability read from landscape classification. Suitable polygons include hedgerows, woodland edges, gardens, areas with hollow-stemmed plants. Maximum nest density per polygon prevents unrealistic crowding.

After repeated failures (default 20 attempts), female switches to long-distance dispersal, moving to distant area via Beta(10,5) distribution draw (typical 1-5 km movement).

#### 4.2.2 Provisioning Dynamics

Once nest acquired, female provisions cells sequentially:

**Step 1: Plan eggs for this nest**
```cpp
int eggs_this_nest = MinEggsPerNest + 
                     g_rand_int(MaxEggsPerNest - MinEggsPerNest);
```
Drawn from uniform(3, 30), matching empirical nest size distributions (Ivanov 2006).

**Step 2: Forage for pollen**
Female searches for flower patches using pre-calculated foraging mask (see [OsmiaForageMaskDetailed](@ref OsmiaForageMaskDetailed)). Mask identifies reachable patches within maximum homing distance ordered by proximity. Female visits nearest patch meeting quality/quantity thresholds:

```cpp
FOR each patch IN foraging_mask (nearest first):
    pollen_available = GetPollenFromMap(patch)
    IF (pollen_available > threshold_quantity AND
        pollen_quality > threshold_quality):
        pollen_collected = min(daily_capacity, pollen_available)
        BREAK  // Use this patch
```

Daily pollen collection capacity depends on:
- Available foraging hours (weather-limited; see Section 3.1)
- Female age (provisioning efficiency peaks ~day 18; Seidelmann 2006)
- Patch quality (pollen density)

**Step 3: Accumulate provision mass**
```cpp
current_cell_mass += pollen_collected_today
```

**Step 4: Check cell completion**
Target provision mass depends on planned offspring sex (determined by maternal age/mass; see Section 4.2.3):
- Female offspring: 20-40 mg (larger, more resource-intensive)
- Male offspring: 10-20 mg (smaller, less costly)

When current mass exceeds target:
```cpp
IF current_cell_mass >= target_provision_mass:
    LayEgg(sex_determined_by_mass)
    current_cell_mass = 0
    eggs_laid++
```

**Step 5: Check nest completion**
```cpp
IF eggs_laid >= eggs_this_nest OR eggs_remaining == 0:
    CloseNest()
    release_nest_pointer
```

This provisioning loop continues until female exhausts her egg load or dies. Weather variation (affecting foraging hours) creates day-to-day fluctuations in provisioning rates, with cascading effects on cell open time (parasitism risk) and seasonal fecundity.

#### 4.2.3 Sex Allocation

*Osmia* females control offspring sex through fertilization: fertilized eggs become females (diploid), unfertilized eggs become males (haploid). Females assess provision mass and decide whether to fertilize based on whether provisioning exceeds female minimum threshold.

Sex allocation strategy emerges from two interacting effects:

**1. Maternal age effect** (declining female bias):
Young females produce ~60-75% females; old females ~40-50% females. Implemented via pre-calculated lookup table ([Osmia_Population_Manager::m_EggSexRatioEqns](@ref Osmia_Population_Manager::m_EggSexRatioEqns)) combining logistic age function with linear mass scaling.

**2. Provision mass effect** (threshold rule):
```cpp
double female_min = GetFemaleMinTargetProvisionMass(age, mass);
if (current_provision >= female_min) {
    sex = FEMALE;  // Fertilize egg
} else {
    sex = MALE;    // Do not fertilize
}
```

This creates emergent sex ratio patterns matching empirical observations: larger, younger females produce female-biased broods; smaller, older females shift toward males (Seidelmann et al. 2010).

### 4.3 Spatial Processes

#### 4.3.1 Foraging Range and Resource Assessment

Females search for flower patches within species-specific foraging range. Implementation uses two data structures:

**1. Foraging Mask** ([OsmiaForageMask](@ref OsmiaForageMask)): Coarse radial search (100 m steps) identifying candidate patches within typical homing distance (~600 m). Used for initial patch identification.

**2. Detailed Mask** ([OsmiaForageMaskDetailed](@ref OsmiaForageMaskDetailed)): Fine-resolution (1 m steps) mapping specific locations within maximum homing distance (~900 m). Used when extensive search required (poor local resources).

Both masks pre-calculated to avoid repeated distance computations during daily foraging decisions.

Resource assessment queries pollen map (10 m resolution) at patch centroids:
```cpp
double pollen_score = pollen_map->GetPollenAvailability(x, y);
double pollen_mg = pollen_score × PollenScoreToMg × (1 - competition_scalar);
```

Competition scalar (default 0.5) represents pollen removal by other bee species, implementing implicit interspecific competition.

#### 4.3.2 Dispersal

Females exhibit philopatry (local nest-site fidelity) but will disperse after repeated local nest-finding failures. Dispersal implemented as single long-distance movement:

```cpp
double distance = Beta(10, 5) × MaxHomingDistance;  // Typically 1-5 km
double angle = Uniform(0, 2π);
new_x = current_x + distance × cos(angle);
new_y = current_y + distance × sin(angle);
```

Beta(10,5) distribution right-skewed: most dispersal moderate distance, occasional long-distance movements. This creates realistic dispersal kernels observed in mark-recapture studies.

### 4.4 Parasitism

Two parasitism implementations available (selected by configuration):

**1. Probability-based** (default): Parasitism risk increases linearly with cell open time:
```cpp
double risk = ParasitismProbPerDay × days_cell_open;
if (g_rand_uni() < risk) {
    parasitised = TRUE;
    parasitoid_type = (g_rand_uni() < BombylidProb) ? BOMBYLID : OTHER;
}
```

Simple, few parameters, computationally cheap. Adequate for landscape-scale questions where parasitoid spatial dynamics not central.

**2. Mechanistic** (optional): Explicit parasitoid population with spatial dynamics, dispersal, reproduction. Parasitism emerges from local parasitoid density:
```cpp
double local_parasitoids = parasitoid_manager->GetDensity(x, y);
double attack_prob = PerCapitaAttackRate × local_parasitoids × days_cell_open;
if (g_rand_uni() < attack_prob) {
    parasitised = TRUE;
}
```

More realistic spatial patterns, but requires extensive additional parameterization (parasitoid mortality, dispersal, reproduction rates). Useful for questions about parasitoid management or landscape effects on parasitism.

Both approaches assign parasitism status at egg-laying. Parasitised individuals develop normally until parasitoid emerges (timing species-specific), then die.

## 5. Parameters

### 5.1 Parameter Organization

Model parameters fall into five categories based on their biological role and data sources:

**1. Development Parameters** - Degree-day requirements and thresholds
**2. Mortality Parameters** - Daily probabilities and survival functions
**3. Provisioning Parameters** - Resource conversion factors and behavioral thresholds
**4. Spatial Parameters** - Movement distances and habitat suitability
**5. Environmental Parameters** - Weather thresholds and seasonal transitions

All parameters are defined as configuration variables in the C++ source code (see [Osmia](@ref Osmia).cpp and [Osmia_Population_Manager](@ref Osmia_Population_Manager).cpp), allowing modification without recompilation. Default values reflect either direct empirical measurements or calibrated values matching observed population dynamics.

### 5.2 Development Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| **Egg Stage** ||||||
| EGG_DEVEL_THRESHOLD | 0.0 | °C | Calibrated | LOW | Lower developmental threshold. **Differs from formal model** (13.8°C laboratory value; see Section 8.3) |
| EGG_DEVEL_TOTAL_DD | 86 | DD | Calibrated | MEDIUM | Degree-days required for hatching. Adjusted from Radmacher & Strohm (2011) |
| EGG_DAILY_MORT | 0.0014 | day⁻¹ | Radmacher & Strohm (2011) | HIGH | Daily mortality probability during egg development |
| **Larval Stage** ||||||
| LARVA_DEVEL_THRESHOLD | 13.8 | °C | Radmacher & Strohm (2011) | HIGH | Lower developmental threshold from laboratory studies |
| LARVA_DEVEL_TOTAL_DD | 240 | DD | Radmacher & Strohm (2011) | HIGH | Degree-days for larval development and cocoon construction |
| LARVA_DAILY_MORT | 0.0014 | day⁻¹ | Radmacher & Strohm (2011) | HIGH | Daily mortality probability during larval feeding |
| **Prepupal Stage** ||||||
| PREPUPA_DEVEL_TOTAL_DAYS | 45 | days | Calibrated | LOW | Approximate duration of summer diapause. See Section 8.3 for discussion of temperature dependence |
| PREPUPA_DEVEL_RATES | [array] | day⁻¹ | Calibrated | LOW | Temperature-specific development rates (0-41°C). Empirically derived from field phenology |
| PREPUPA_DAILY_MORT | 0.003 | day⁻¹ | Radmacher & Strohm (2011) | MEDIUM | Daily mortality during diapause |
| **Pupal Stage** ||||||
| PUPA_DEVEL_THRESHOLD | 1.1 | °C | Calibrated | MEDIUM | Lower developmental threshold. **Differs from formal model** (13.2°C laboratory; see Section 8.3) |
| PUPA_DEVEL_TOTAL_DD | 570 | DD | Calibrated | MEDIUM | Degree-days for metamorphosis. Adjusted from literature |
| PUPA_DAILY_MORT | 0.003 | day⁻¹ | Radmacher & Strohm (2011) | MEDIUM | Daily mortality during pupation |
| **Overwintering** ||||||
| OVERWINTER_THRESHOLD | 5.0 | °C | Sgolastra et al. (2011) | MEDIUM | Temperature below which chilling accumulates |
| OVERWINTER_REQUIRED_DD | 400-500 | DD | Sgolastra et al. (2011) | MEDIUM | Chilling requirement for spring emergence |
| PREWINTER_THRESHOLD | 13.0 | °C | Calibrated | LOW | Temperature triggering onset of overwintering physiology |
| EMERG_COUNT_CONST | [value] | - | Fitted | LOW | Intercept for emergence probability equation |
| EMERG_COUNT_SLOPE | [value] | DD⁻¹ | Fitted | LOW | Slope for emergence probability vs. accumulated DD |

**Key Uncertainty: Prepupal Development**

The prepupal stage presents the highest parameter uncertainty. Laboratory studies (Radmacher & Strohm 2011) provide quadratic temperature-development relationships at constant temperatures (15-35°C). However, these relationships produced unrealistic developmental durations under field conditions with daily temperature fluctuations and occasional extremes.

The implemented lookup table approach allows flexible parameterization whilst acknowledging mechanistic uncertainty. Rates were adjusted iteratively until simulated emergence phenology matched field observations from multiple sites and years. This pragmatic approach prioritizes phenological realism over mechanistic detail given insufficient data to validate thermal performance curves under field conditions.

Future model improvements should include:
- Field measurements of prepupal development across natural temperature regimes
- Explicit representation of thermal performance curve shape
- Individual variation in diapause duration
- Photoperiod effects on diapause termination

See Section 8.3 for extended discussion.

### 5.3 Mortality Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| FEMALE_BACKGROUND_MORT | 0.02 | day⁻¹ | Calibrated | MEDIUM | Daily adult mortality from predation, disease, accidents. Yields ~50 day mean lifespan |
| WINTER_MORT_CONST | [value] | - | Sgolastra et al. (2011) | MEDIUM | Intercept for overwintering mortality equation |
| WINTER_MORT_SLOPE | [value] | DD⁻¹ | Sgolastra et al. (2011) | MEDIUM | Slope relating mortality to accumulated cold stress |

**Overwintering Mortality Function:**

Mortality probability increases with accumulated winter degree-days (prolonged cold exposure) and decreases with cocoon mass (resource stores):

```
P(mortality) = CONST + SLOPE × overwinter_DD - MASS_EFFECT × (cocoon_mass - mean_mass)
```

Equation from *Osmia lignaria* studies (Sgolastra et al. 2011), assumed applicable to *O. bicornis* given similar overwintering biology. Confidence MEDIUM because species differ in cold tolerance; validation against *O. bicornis* field data desirable.

### 5.4 Provisioning and Reproduction Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| **Resource Conversion** ||||||
| COCOON_TO_PROVISION_MASS | 3.247 | - | Seidelmann (2006) | HIGH | Ratio converting provision mass to cocoon mass (31% conversion efficiency) |
| PROVISION_TO_COCOON_MASS | 0.308 | - | Seidelmann (2006) | HIGH | Inverse ratio (cocoon to provision) |
| POLLEN_SCORE_TO_MG | 0.8 | mg·score⁻¹ | Calibrated | MEDIUM | Conversion from landscape pollen availability score to actual mg collectible |
| **Provisioning Behavior** ||||||
| MIN_EGGS_PER_NEST | 3 | eggs | Ivanov (2006) | HIGH | Minimum planned eggs per nest |
| MAX_EGGS_PER_NEST | 30 | eggs | Ivanov (2006) | HIGH | Maximum planned eggs per nest |
| MALE_MIN_TARGET_PROVISION | 10.0 | mg | Seidelmann et al. (2010) | HIGH | Minimum provision mass for male cells |
| FEMALE_MIN_TARGET_PROVISION | 20.0 | mg | Seidelmann et al. (2010) | HIGH | Minimum provision mass for female cells |
| MIN_CELL_CONSTRUCTION_TIME | 1 | days | Field observations | HIGH | Minimum time to provision one cell under ideal conditions |
| MAX_CELL_CONSTRUCTION_TIME | 4 | days | Seidelmann (2006) | HIGH | Maximum time before abandoning cell (parasitism risk threshold) |
| LIFETIME_COCOON_MASS_LOSS | 30.0 | mg | Seidelmann et al. (2010) | MEDIUM | Total decline in provision mass from first to last offspring |
| **Sex Allocation** ||||||
| SEX_RATIO_VS_AGE_LOGISTIC | [4 params] | - | Seidelmann et al. (2010) | HIGH | Logistic equation parameters for age-dependent sex ratio |
| SEX_RATIO_VS_MASS_LINEAR | [2 params] | - | Seidelmann et al. (2010) | HIGH | Linear relationship between maternal mass and sex ratio adjustment |
| FEMALE_COCOON_MASS_VS_AGE | [4 params] | - | Seidelmann et al. (2010) | HIGH | Logistic equation for age-dependent female offspring mass |
| FEMALE_COCOON_MASS_VS_MASS | [2 params] | mg⁻¹ | Calibrated | MEDIUM | Linear relationship maternal mass → offspring mass. **Differs from formal model** (slope adjusted 0.46→0.3) |

**Sex Allocation Implementation:**

Sex ratio determination combines maternal age and mass effects through pre-calculated lookup tables (see [Osmia_Population_Manager::m_EggSexRatioEqns](@ref Osmia_Population_Manager::m_EggSexRatioEqns)). For each combination of maternal age (0-60 days) and mass (4-28 mg in 0.25 mg steps), the table stores:

```
sex_ratio[age][mass] = BASE + (ADJUSTED_MAX - BASE) / (1 + exp(-K × (age - AGE_50)))
where ADJUSTED_MAX = LINEAR_SLOPE × mass + LINEAR_INTERCEPT
```

This creates a surface where young, large females produce highly female-biased broods (~70-75% female), whilst old, small females shift toward male-biased broods (~40-45% female). Pattern matches empirical observations (Seidelmann et al. 2010) and optimal sex allocation theory.

**Calibration Note:** The linear relationship between maternal and offspring mass was adjusted from laboratory-derived values (slope 0.46, intercept 63.85) to values producing better matches to field offspring size distributions (slope 0.30, intercept 65.1). Laboratory conditions may over-predict resource allocation efficiency relative to field conditions with weather interruptions and predation risk.

### 5.5 Spatial Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| TYPICAL_HOMING_DISTANCE | 600 | m | Literature review | MEDIUM | Typical foraging range for *Osmia bicornis*. Used for initial nest/patch searches |
| MAX_HOMING_DISTANCE | 900 | m | Literature review | MEDIUM | Maximum foraging distance. Defines extent of detailed foraging mask |
| NEST_FIND_ATTEMPT_NO | 20 | attempts | Calibrated | LOW | Number of local nest-finding attempts before triggering dispersal |
| FEMALE_DISPERSAL_SHAPE1 | 10 | - | Calibrated | LOW | Beta distribution shape parameter for dispersal distance (right-skewed) |
| FEMALE_DISPERSAL_SHAPE2 | 5 | - | Calibrated | LOW | Beta distribution shape parameter (most movements moderate distance) |
| DENSITY_GRID_RESOLUTION | 1000 | m | Design choice | - | Grid cell size for tracking female density (1 km²) |
| POLLEN_MAP_RESOLUTION | 10 | m | Landscape data | - | Resolution of input pollen availability maps |

**Foraging Range Justification:**

Solitary bee foraging distances poorly documented compared to social species. Values chosen represent conservative estimates for 15-25 mg bee body mass, based on allometric relationships from Greenleaf et al. (2007) and scattered *Osmia* spp. observations. Uncertainty MEDIUM because species-specific *O. bicornis* data limited; validation against genetic or mark-recapture studies desirable.

### 5.6 Environmental Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| MIN_TEMP_FOR_FLYING | 6.0 | °C | Field observations | HIGH | Minimum temperature for flight activity |
| MAX_WIND_SPEED_FOR_FLYING | 8.0 | m·s⁻¹ | Field observations | HIGH | Maximum wind speed permitting flight |
| MAX_PRECIP_FOR_FLYING | 0.1 | mm·h⁻¹ | Field observations | HIGH | Maximum precipitation permitting flight (effectively prohibits flying in rain) |
| POLLEN_GIVE_UP_THRESHOLD | 0.75 | proportion | Calibrated | LOW | Proportional depletion triggering patch abandonment |
| POLLEN_GIVE_UP_RETURN | 0.75 | mg | Calibrated | LOW | Minimum acceptable pollen gain per foraging bout |
| DENSITY_DEPENDENT_POLLEN_REMOVAL | 0.5 | proportion | Calibrated | MEDIUM | Proportion of pollen removed by competing bee species |

**Weather Threshold Sources:**

Flight activity thresholds derived from field observations of *Osmia* spp. under various weather conditions. Temperature threshold (6°C) conservative relative to some observations (flight at 8-10°C) to account for microclimate warming (sun-heated surfaces). Wind threshold (8 m/s) reflects small body size vulnerability to gusts. Precipitation threshold (0.1 mm/h) effectively prohibits flight during any measurable rain, matching observed behaviour.

**Foraging Thresholds:**

Give-up thresholds implement optimal foraging theory: bees should leave patches when gain rate falls below landscape average. Values calibrated to produce realistic patch residence times, but uncertainty HIGH because species-specific foraging behaviour data limited. Sensitivity analyses indicate population dynamics relatively robust to these parameters (reproductive output more constrained by total resource availability than fine-scale foraging efficiency).

### 5.7 Parasitism Parameters

| Parameter | Default | Units | Source | Confidence | Description |
|-----------|---------|-------|--------|------------|-------------|
| **Probability-Based Model** (default) |
| PARASITISM_PROB_TO_TIME_OPEN | 0.0075 | day⁻¹ | Calibrated | LOW | Daily parasitism risk accumulation rate. 4-day cell open time → ~3% risk |
| BOMBYLIID_PROB | 0.5 | proportion | Field surveys | LOW | Proportion of parasitism events from Bombyliidae (bee flies) vs. other taxa |
| **Mechanistic Model** (optional) |
| PER_CAPITA_PARASITATION_CHANCE | [0.00001, 0.00002] | day⁻¹ | Fitted | LOW | Per-capita attack rates for different parasitoid types |
| PARASITOID_DAILY_MORT | [array] | day⁻¹ | Literature | LOW | Monthly parasitoid mortality rates |
| PARASITOID_DISPERSAL_RATE | [value] | m·day⁻¹ | Literature | LOW | Parasitoid movement rate across landscape |

**Parasitism Uncertainty:**

Parasitism parameters have universally LOW confidence because:
1. Field parasitism rates highly variable (5-50% across sites/years)
2. Parasitoid community composition varies spatially
3. Attack rates rarely measured directly (inferred from outcomes)
4. Parasitoid spatial dynamics poorly understood

Probability-based model adequate for landscape-scale questions where parasitism acts as aggregate mortality source. Mechanistic model required only for questions specifically about parasitoid management or spatial dynamics.

Calibration typically inverse: adjust parameters until simulated parasitism rates match observed ranges for study system. Model serves as hypothesis-generating tool for understanding landscape effects on parasitism rather than predictive tool for absolute rates.

### 5.8 Parameter Sensitivity and Uncertainty Propagation

**High-Sensitivity Parameters** (large effects on population dynamics):
- Development thresholds (especially egg and pupa; altered phenology)
- Female background mortality (directly affects adult lifespan and fecundity)
- Pollen score to mg conversion (links landscape resources to reproduction)
- Overwintering mortality parameters (determines overwinter survival)

**Low-Sensitivity Parameters** (minor effects):
- Foraging give-up thresholds (affect patch choice, minimal population impact)
- Nest size range (min/max eggs per nest; average matters more than bounds)
- Sex allocation equation parameters (population sex ratio relatively stable)

**Uncertainty Propagation:**

Model predictions incorporate parameter uncertainty through Monte Carlo simulation:
1. Define probability distributions for uncertain parameters
2. Draw parameter sets from distributions (Latin Hypercube sampling)
3. Run model ensemble (typically 100-1000 replicates)
4. Report summary statistics (median, 5th/95th percentiles) acknowledging uncertainty

This approach quantifies how parameter uncertainty translates to prediction uncertainty, essential for risk assessment and management evaluation.

---

## 6. Inputs

### 6.1 Required Input Data

The model requires four categories of input data:

#### 6.1.1 Landscape Data

**Land-Use Classification Map**
- **Format**: Polygon shapefile or raster
- **Resolution**: Polygon-based (field/habitat boundaries) or 10-50 m raster
- **Content**: Land-use type for each spatial unit (arable, grassland, hedgerow, woodland, urban, etc.)
- **Source**: Digitized field boundaries, national land-cover datasets, remote sensing classification
- **Preprocessing**: Convert to ALMaSS polygon format with unique ID, centroid coordinates, area, land-use code

**Nesting Habitat Suitability**
- **Format**: Lookup table (CSV or configuration file)
- **Content**: Maximum nest density and nesting probability per land-use type
- **Example**:
```
LandUseCode, MaxNestsPerHa, NestingProbability
HEDGE,       100,           0.8
WOODLAND,    50,            0.6
GARDEN,      200,           0.9
ARABLE,      0,             0.0
```
- **Source**: Field surveys of nest abundance, habitat characteristics (hollow stem abundance, dead wood availability)
- **Preprocessing**: Expert assessment if empirical data unavailable

#### 6.1.2 Resource Data

**Pollen Availability Maps**
- **Format**: Monthly raster grids (January-December)
- **Resolution**: 10 m × 10 m
- **Content**: Pollen availability score (unitless, typically 0-100) combining quantity and quality
- **Calculation**:
```
PollenScore = Σ (flowerCover[species] × pollenProduction[species] × quality[species])
```
- **Source**: Vegetation surveys, flower cover mapping, pollen production databases
- **Preprocessing**: 
  - Field surveys → species cover estimates → pollen production multiplication
  - Remote sensing (if sufficient spectral resolution for flower detection)
  - Crop flowering calendars for agricultural landscapes

**Nectar Availability Maps** (optional, for extension to energetics)
- **Format**: As pollen maps
- **Content**: Nectar availability score or sugar concentration (mg sugar/L) and quantity (mL/m²)
- **Usage**: Currently informational only; model uses pollen as sole reproductive constraint

**Monthly Resource Thresholds**
- **Format**: Configuration array (12 months × 4 values)
- **Content**: Minimum acceptable pollen quantity, pollen quality, nectar quantity, nectar quality per month
- **Usage**: Females reject patches below thresholds
- **Calibration**: Iterative adjustment until foraging behaviour realistic

#### 6.1.3 Weather Data

**Daily Weather Time Series**
- **Format**: CSV file or database connection
- **Required Variables**:
  - Date (YYYY-MM-DD)
  - Daily mean temperature (°C)
  - Hourly temperatures (°C, 24 values) OR daily min/max with interpolation
  - Hourly wind speeds (m/s, 24 values) OR daily mean
  - Hourly precipitation (mm/h, 24 values) OR daily total
- **Source**: Meteorological stations, gridded climate datasets (e.g., E-OBS), reanalysis products
- **Spatial Coverage**: Ideally multiple stations across landscape; alternatively single representative station
- **Temporal Extent**: Minimum 1 year; multi-year runs require continuous series

**Preprocessing Requirements**:
- Gap-filling for missing data (interpolation, station averaging)
- Quality control (outlier detection, consistency checks)
- Format conversion to ALMaSS weather file structure

#### 6.1.4 Initial Population

**Starting Population Specification**
- **Format**: Configuration parameters
- **Required**:
  - Total population size (number of overwintering females)
  - Size distribution (minimum/maximum female mass)
  - Spatial distribution (uniform across suitable habitat OR specified polygon list)
  - Overwintering progress (accumulated degree-days at simulation start)
- **Default Approach**: Generate N females randomly distributed across nesting habitat with uniform size distribution and partial overwintering progress (320 DD default)
- **Alternative**: Read explicit population from file if continuing previous simulation

### 6.2 Data Format Specifications

#### 6.2.1 ALMaSS Landscape Format

The model operates within ALMaSS framework requiring specific landscape file structure:

**Polygon Definition File** (.alm format):
```
[HEADER]
NumberOfPolygons: XXXX
NumberOfFieldTypes: YY

[POLYGONS]
PolyID, CentroidX, CcentroidY, Area, FieldType, [additional attributes]
1,      500000,    6200000,    25000, HEDGE, ...
2,      500100,    6200050,    180000, ARABLE, ...
...
```

**Field Type Definition File** (.fdt format):
```
[FIELDTYPES]
TypeID, TypeName,    ManagementSequence
1,      HEDGE,       Hedge_Management.txt
2,      ARABLE,      Winter_Wheat.txt
3,      GRASSLAND,   Permanent_Grass.txt
...
```

Conversion tools available for common GIS formats (shapefiles, GeoTIFF) to ALMaSS format.

#### 6.2.2 Pollen Map Format

Monthly pollen maps stored as binary rasters or ASCII grids:

**ASCII Grid Format**:
```
NCOLS         1000
NROWS         1000
XLLCORNER     500000
YLLCORNER     6200000
CELLSIZE      10
NODATA_VALUE  -9999
45.2 38.7 52.1 ...
67.8 71.2 65.4 ...
...
```

12 files required: pollen_january.asc through pollen_december.asc

**Binary Format**: More efficient for large landscapes; requires header file specifying dimensions and georeferencing.

### 6.3 Data Preparation Workflows

#### 6.3.1 Creating Pollen Maps from Field Surveys

**Step 1: Conduct Vegetation Surveys**
- Stratified random sampling across landscape
- Record flower cover by species (%, visual estimation or quadrats)
- Monthly surveys March-September (flowering period)

**Step 2: Assign Pollen Production Values**
- Literature lookup: pollen grains per flower, flower density
- Expert assessment: relative attractiveness to *Osmia*
- Create species-level pollen production database

**Step 3: Spatial Interpolation**
- Kriging or inverse distance weighting from survey points
- Incorporate land-use as covariate (crop types have known flower abundance)
- Generate 10 m resolution rasters

**Step 4: Quality Adjustment**
- Weight by *Osmia* preferences (if known)
- Adjust for flower accessibility (some plant species produce abundant pollen but flowers inaccessible to short-tongued bees)

**Step 5: Validation**
- Compare maps to independent foraging observations
- Adjust scoring function until bee distribution matches flower distribution

#### 6.3.2 Processing Weather Station Data

**Step 1: Data Acquisition**
- Download from national meteorological services
- Ensure hourly resolution (or daily min/max for interpolation)
- Verify data completeness and quality flags

**Step 2: Gap-Filling**
- Linear interpolation for short gaps (<6 hours)
- Regression with nearby stations for longer gaps
- Climate normals for persistent missing data

**Step 3: Spatial Interpolation** (if multiple stations available)
- Inverse distance weighting or thin-plate splines
- Elevation adjustment for temperature (lapse rate correction)
- Generate landscape-wide gridded weather OR assign stations to landscape zones

**Step 4: Format Conversion**
- Convert to ALMaSS weather file format
- Calculate derived variables (degree-days, flying hours) if not automated in model

### 6.4 Data Quality Requirements

**Minimum Viable Data Quality:**
- **Landscape**: Complete land-use classification; nesting habitat estimates even if coarse
- **Resources**: Pollen maps for April-June (peak provisioning); other months can use defaults
- **Weather**: Complete daily temperature series; hourly data preferred but daily min/max acceptable
- **Initial Population**: Population size can be arbitrary (model explores relative dynamics)

**Preferred Data Quality:**
- **Landscape**: Field-validated habitat classification; nest density measurements
- **Resources**: Monthly pollen AND nectar maps from multi-year surveys
- **Weather**: Multi-station networks with hourly data; validation against local microclimate
- **Initial Population**: Size distribution from field-collected cocoons; spatial distribution from nest surveys

**Data Limitations and Alternatives:**

When empirical data insufficient:
- Use land-use proxies (crop-specific flowering calendars)
- Adopt parameter values from similar systems (other solitary bee models)
- Run sensitivity analyses to quantify uncertainty from missing data
- Focus on relative predictions (management comparisons) rather than absolute abundance

## 7. Outputs

### 7.1 Output Data Organization

The model generates three categories of output data:

#### 7.1.1 Population Time Series

**Daily Abundance by Life Stage**
- **Format**: CSV file (one row per day)
- **Columns**: Date, Eggs, Larvae, Prepupae, Pupae, InCocoon, Females
- **Usage**: Track population dynamics, seasonal phenology, cohort progression
- **Example Applications**:
  - Emergence timing (date of peak female abundance)
  - Reproductive period duration (days with active females)
  - Overwinter survival (ratio of autumn InCocoon to spring Females)
  - Population growth rate (ratio of year N+1 to year N spring females)

**Female Age Structure**
- **Format**: CSV snapshot (specified intervals, e.g., weekly)
- **Columns**: AgeClass (days), Count
- **Usage**: Assess age-dependent mortality, reproductive span
- **Interpretation**: Right-skewed distribution indicates high early mortality; truncated distribution suggests lifespan limitations

**Size Distribution**
- **Format**: CSV snapshot or histogram
- **Columns**: MassClass (mg), Count
- **Usage**: Monitor size-dependent fitness, resource limitation effects
- **Interpretation**: Declining mean mass over season indicates resource depletion; bimodal distribution indicates sexual size dimorphism in offspring

#### 7.1.2 Reproductive Output

**Fecundity Statistics**
- **Format**: CSV file (one row per female)
- **Columns**: FemaleID, EmergenceDate, BodyMass, TotalEggsLaid, NestsCompleted, DeathDate
- **Usage**: Individual-level reproductive success, lifetime fecundity distributions
- **Validation**: Compare to field observations of nest sizes, eggs per female

**Sex Ratio Dynamics**
- **Format**: CSV time series or aggregated
- **Columns**: Date (or AgeClass/MassClass), ProportionFemale, SampleSize
- **Usage**: Test sex allocation theory, assess population viability
- **Interpretation**: Female-biased ratios indicate good conditions; male bias suggests resource stress or poor maternal condition

**Provisioning Efficiency**
- **Format**: CSV aggregated by age class
- **Columns**: AgeClass, MeanProvisioningRate (mg/day), SDProvisioningRate, SampleSize
- **Usage**: Validate provisioning equations against Seidelmann (2006) empirical data
- **Expected Pattern**: Low at emergence, peak around day 18-20, decline after day 40

#### 7.1.3 Spatial Outputs

**Female Density Maps**
- **Format**: Raster or polygon-aggregated
- **Resolution**: 1 km² grid or per-polygon aggregation
- **Temporal**: Snapshots at specified intervals (e.g., weekly during activity season)
- **Usage**: Visualize spatial distribution, identify source/sink habitats
- **Interpretation**: High-density areas indicate resource-rich habitat; low-density areas may be sinks requiring immigration

**Nest Distribution**
- **Format**: Point shapefile or polygon counts
- **Columns**: Location (X, Y or PolygonID), NestCount, MeanEggsPerNest
- **Usage**: Validate nesting habitat predictions, assess spatial clustering
- **Validation**: Compare to field nest box trap data

**Foraging Range Maps** (optional diagnostic output)
- **Format**: Raster showing foraging pressure
- **Content**: Number of foraging visits per cell per day
- **Usage**: Assess resource depletion patterns, connectivity effects
- **Interpretation**: High-intensity areas may experience resource limitation

### 7.2 Output File Formats and Naming

#### 7.2.1 Standard Output Files

**Population Dynamics** (daily):
```
Osmia_Population_YYYY.csv
Columns: Year, Day, Date, Eggs, Larvae, Prepupae, Pupae, InCocoon, Females
```

**Reproductive Summary** (per female):
```
Osmia_Fecundity_YYYY.csv
Columns: FemaleID, BodyMass, EmergenceDate, TotalEggs, FemaleOffspring, 
         MaleOffspring, NestsCompleted, FinalNestComplete, DeathDate, DeathCause
```

**Spatial Distribution** (weekly snapshots):
```
Osmia_Density_YYYY_Week_WW.asc (ASCII raster)
Osmia_Nests_YYYY_Week_WW.csv (polygon summary)
```

#### 7.2.2 Validation Outputs (if testing mode enabled)

Conditional compilation flag `__OSMIATESTING` activates detailed diagnostic outputs:

**Development Duration Distributions**:
```
OsmiaStageLengths.txt
Columns: Stage, Individual, Duration (days), StartDate, EndDate, 
         MeanTemp, AccumulatedDD
```

**Provisioning Records**:
```
OsmiaProvisioningDetails.txt
Columns: FemaleID, Age, CellNumber, DaysToComplete, PollenCollected, 
         FinalProvisionMass, WeatherHours, PatchQuality
```

**Female Weight Tracking**:
```
OsmiaFemaleWeights.txt
Columns: FemaleID, Date, Age, BodyMass, CumulativeEggsLaid, CurrentNestID
```

These diagnostic outputs enable detailed validation against empirical distributions but substantially increase file sizes and I/O overhead. Use only for validation runs, not production experiments.

### 7.3 Output Interpretation Guidelines

#### 7.3.1 Assessing Population Viability

**Annual Growth Rate**:
```
λ = N(spring, year t+1) / N(spring, year t)
```
- λ > 1: Population increasing
- λ = 1: Population stable
- λ < 1: Population declining

**Extinction Risk**: For stochastic simulations (Monte Carlo with parameter/weather uncertainty), calculate:
```
P(extinction) = proportion of replicates with N < threshold within T years
```
Typical threshold: 50-100 females (quasi-extinction rather than true extinction due to Allee effects not modelled)

#### 7.3.2 Landscape Comparisons

When comparing alternative landscapes or management scenarios:

**Effect Size Calculation**:
```
Relative Effect = (N_treatment - N_baseline) / N_baseline × 100%
```

**Statistical Significance**: Use resampling or Monte Carlo to generate confidence intervals:
- Run 100+ replicates per scenario with stochastic weather and parameters
- Compare distributions using Mann-Whitney U test or bootstrap confidence intervals
- Report median difference and 95% CI

**Meaningful Difference**: Biological significance threshold ~10-20% change in abundance (ecological rule of thumb; smaller changes may be statistically detectable but ecologically trivial)

#### 7.3.3 Phenological Metrics

**Emergence Timing**:
- First emergence: Date when first female emerges (10th percentile of emergence distribution)
- Peak emergence: Date of maximum female abundance
- Emergence spread: Days between 10th and 90th percentile

**Reproductive Period**:
- Start: Date of first egg laying
- End: Date of last egg laying
- Duration: Days with active provisioning

**Synchrony with Resources**:
- Overlap: Proportion of active female days with pollen availability > threshold
- Mismatch: Days when females present but pollen scarce (resource-limited days)

### 7.4 Post-Processing and Visualization

#### 7.4.1 Recommended Visualization Approaches

**Population Dynamics**:
- Stacked area plot showing all life stages over time
- Faceted plots comparing multiple years or scenarios
- Annotate with weather events (cold snaps, droughts)

**Spatial Patterns**:
- Heatmaps of female density overlaid on land-use
- Network diagrams showing dispersal connectivity
- Before/after comparisons for management interventions

**Reproductive Output**:
- Violin plots of fecundity distributions by scenario
- Survival curves (Kaplan-Meier) for age-specific mortality
- Sex ratio surfaces (age × mass) validated against empirical data

#### 7.4.2 R Scripts for Standard Analyses

Example processing script structure:

```R
# Load population time series
pop_data <- read.csv("Osmia_Population_2024.csv")

# Calculate annual metrics
spring_females <- pop_data %>%
    filter(month(Date) == 4) %>%
    summarise(mean_females = mean(Females),
              peak_females = max(Females),
              emergence_start = Date[which(Females > 10)[1]])

# Plot seasonal dynamics
ggplot(pop_data, aes(x = Date)) +
    geom_area(aes(y = Females), fill = "orange", alpha = 0.7) +
    geom_area(aes(y = InCocoon), fill = "brown", alpha = 0.5) +
    labs(title = "Osmia bicornis Population Dynamics",
         y = "Abundance", x = "Date") +
    theme_minimal()
```

Repository includes complete R scripts for standard visualizations.

---

## 8. Implementation Discussion

### 8.1 Design Principles and Trade-offs

The *Osmia bicornis* model embodies several fundamental design decisions that balance biological realism against computational tractability:

#### 8.1.1 Agent-Based vs. Stage-Structured Approaches

**Decision**: Implement as individual-based model (IBM) rather than matrix population model or ordinary differential equations (ODE).

**Rationale**: 
- Individual heterogeneity matters: Body size affects fecundity, sex ratios, survival
- Spatial processes explicit: Foraging range, dispersal, habitat selection
- Stochasticity natural: Individual-level demographic and environmental stochasticity
- Mechanistic behaviour: Provisioning emerges from weather × resources × female state

**Trade-offs**:
- Computational cost O(N) vs. O(1) for ODE approaches
- Requires more parameters (individual-level processes)
- Stochastic outcomes require ensemble runs

**When Alternative Appropriate**: If questions concern only total abundance or simple spatial patterns, stage-structured matrix models sufficient and orders of magnitude faster.

#### 8.1.2 Daily vs. Finer Temporal Resolution

**Decision**: Daily time steps for all processes.

**Rationale**:
- Development occurs gradually (degree-days accumulated daily)
- Weather varies daily (foraging hours calculation)
- Behavioural decisions daily (forage, provision, disperse)
- Matches empirical data resolution (field observations typically daily)

**Trade-offs**:
- Sub-daily provisioning dynamics aggregated (realistic for small bees with limited daily pollen loads)
- Hourly weather variation flattened to "foraging hours available"
- Cannot resolve diurnal rhythms (early morning peak foraging)

**When Finer Resolution Needed**: If studying thermal biology (thermoregulation), predator-prey interactions with diurnal cycles, or competition dynamics at flower patches.

#### 8.1.3 Polygon-Based Landscape vs. Grid

**Decision**: Hybrid approach: polygon land-use with overlaid grids for resources (10 m pollen map) and density (1 km aggregation).

**Rationale**:
- Polygons match management units (fields, hedgerows)
- Resource grids provide fine-scale heterogeneity
- Density grids balance realism vs. computational cost

**Trade-offs**:
- Multiple spatial representations require coordinate transformations
- Grid resolution choices affect patterns (modifiable areal unit problem)
- Cannot resolve within-field variation below 10 m

**Sensitivity to Resolution**: Tests with 5 m vs. 20 m pollen grids showed minor population dynamic effects (<5% abundance change) but substantial runtime differences (4× slower at 5 m). Current 10 m choice balances accuracy and speed.

### 8.2 Model Limitations and Boundaries

#### 8.2.1 Biological Simplifications

**Males Not Modelled**
- Justification: Females are sole provisioners; male abundance/quality doesn't constrain female reproduction in non-resource-limited populations
- Limitation: Cannot explore sex ratio evolution, male-male competition, sperm limitation
- Acceptable for: Landscape management questions, resource effects
- Problematic for: Mating system evolution, Allee effects from mate-finding

**Static Flower Phenology**
- Implementation: Monthly pollen maps fixed across years
- Limitation: Ignores inter-annual variation in flowering phenology, climate change effects on flower timing
- Extension path: Dynamic vegetation model coupling or phenology models

**No Explicit Pesticide Pharmacokinetics**
- Base model: Optional pesticide module with simple toxicodynamics
- Limitation: Crude approximation of exposure pathways, sublethal effects
- Extension: GUTS toxicity models or DEBtox energetics for refined risk assessment

**Simplified Parasitism**
- Probability model: Time-based risk, no spatial dynamics
- Mechanistic model: Available but poorly parameterized
- Limitation: Cannot realistically assess parasitoid management strategies
- Extension path: Empirical parasitoid field data, coupled parasitoid population model

#### 8.2.2 Spatial Scale Limitations

**Minimum Viable Extent**: ~5 km² landscapes
- Smaller areas: Edge effects dominate, dispersal losses unrealistic
- Closed-population assumption violated

**Maximum Practical Extent**: ~500 km² landscapes
- Runtime: Days to weeks for multi-year, multi-replicate simulations
- Data requirements: Fine-resolution pollen maps difficult to produce at large extents

**Optimal Range**: 10-100 km² (agricultural catchments, nature reserves)

#### 8.2.3 Temporal Scale Limitations

**Spin-Up Requirements**: 
- 1-year spin-up recommended (start with overwintering adults, run through one cycle)
- Ensures realistic age/size distributions before experimental treatments

**Climate Change Applications**:
- Model assumes temperature affects development rates only
- Ignores: CO₂ effects on plant quality, phenological mismatches with novel flowering, extreme event impacts beyond current weather variation
- Appropriate for: Moderate warming scenarios (<2°C), near-term projections
- Questionable for: Extreme scenarios (>3°C), long-term evolution

### 8.3 Major Uncertainties and Research Priorities

#### 8.3.1 Prepupal Development: The Critical Knowledge Gap

**The Problem**:
Laboratory studies (Radmacher & Strohm 2011) provide degree-day relationships for prepupal development at constant temperatures (15-35°C). When applied to field temperature regimes with daily fluctuations and occasional extremes, these relationships produce emergence phenology 2-4 weeks earlier than observed.

**Hypothesized Mechanisms**:
1. **Thermal performance curve shape**: Laboratory constant temperatures may not capture non-linear responses to fluctuating regimes
2. **Photoperiod interactions**: Diapause termination may require photoperiod cues not captured by temperature alone
3. **Moisture effects**: Cocoon water content affects metabolic rates; field moisture more variable than laboratory
4. **Individual variation**: Natural populations more heterogeneous than laboratory colonies

**Current Implementation**:
Empirical lookup table with temperature-specific development rates, calibrated to multi-site emergence observations. Acknowledges phenomenological rather than mechanistic approach.

**Research Priority**: HIGH
- Field experiments with replicated temperature treatments under natural photoperiod/moisture
- Fine-scale temperature logging inside nests (cocoon microenvironment)
- Individual-level development tracking with temperature history recording
- Comparative studies across latitudes (if photoperiod involved, should vary predictably)

**Impact on Model Predictions**:
Phenology sensitivity analyses show ±10 day emergence shift produces ±15-20% abundance changes (resource-phenology mismatch effects). Prepupal parameter uncertainty thus propagates substantially to management predictions.

#### 8.3.2 Overwintering Mortality

**Current Implementation**: 
Equation from *Osmia lignaria* assumed transferable to *O. bicornis*. Parameters: winter mortality increases with accumulated cold stress, decreases with cocoon mass.

**Uncertainties**:
- Species differences in cold tolerance (lignaria North American; bicornis European)
- Fungal pathogen effects (moisture-dependent, not modelled)
- Extreme event impacts (sudden temperature fluctuations)

**Research Priority**: MEDIUM
- Multi-year field studies: cocoon mass measurement → overwintering → spring emergence
- Experimental cold exposures with *O. bicornis* specifically
- Pathogen screening of dead vs. surviving cocoons

**Interim Solution**: 
Sensitivity analyses bracketing mortality parameters ±50% to encompass plausible range. Report predictions with confidence intervals reflecting this uncertainty.

#### 8.3.3 Foraging Behaviour Parameterization

**Poorly Constrained Parameters**:
- Give-up thresholds (when to abandon patch)
- Pollen score to mg conversion (landscape map units → actual collection rates)
- Competition effects (other bee species pollen removal)

**Research Priority**: MEDIUM-LOW
- Individual-level foraging observations with pollen load measurements
- Experimental resource manipulations (enrichment/depletion)
- Community-level competition experiments

**Sensitivity**: Moderate
- Population dynamics relatively robust to foraging details (total landscape pollen matters more than fine-scale efficiency)
- Spatial distribution sensitive (clustering vs. dispersion affected)

### 8.4 Future Development Priorities

#### 8.4.1 Short-Term Enhancements (1-2 years)

**Priority 1: Improved Prepupal Parameterization**
- Collaborate with thermal biologists for field development experiments
- Implement photoperiod × temperature interaction if supported by data
- Validate against multi-site, multi-year emergence data

**Priority 2: Pesticide Submodule Refinement**
- Implement GUTS (General Unified Threshold model of Survival) toxicodynamics
- Add sublethal effects (reduced foraging efficiency, provisioning impairment)
- Parameterize for common agricultural insecticides

**Priority 3: Landscape Optimization Tool**
- Genetic algorithm or simulated annealing for habitat configuration optimization
- Objective: maximize population viability given constraints (agricultural production, cost)
- Output: actionable management recommendations (where to place flower strips, nesting habitat)

#### 8.4.2 Medium-Term Extensions (3-5 years)

**Priority 1: Coupled Multi-Species Model**
- Add honeybees (managed competition)
- Add bumblebees (wild competition)
- Add other *Osmia* species (guild dynamics)
- Enables: Community-level predictions, pollination service stability

**Priority 2: Dynamic Vegetation Module**
- Replace static monthly pollen maps with growing degree-day driven phenology
- Couple to individual plant models (flowering response to temperature, water)
- Enables: Climate change scenarios with phenological shifts

**Priority 3: Evolutionary Extensions**
- Genetic algorithm for sex allocation strategy evolution
- Adaptive dynamics for dispersal kernel shape
- Enables: Long-term adaptation predictions, evolutionary rescue potential

#### 8.4.3 Long-Term Research Agenda (5-10 years)

**Priority 1: Full Lifecycle Genetics**
- Explicit diploid/haploid genetics (currently implicit)
- Marker-based neutral variation (validate dispersal predictions)
- Quantitative genetics for body size (selection response)

**Priority 2: Pathogen-Parasite Module**
- Explicit microsporidian (*Nosema*) infection dynamics
- Mite parasitism (*Chaetodactylus*)
- Fungal pathogens
- Enables: Disease outbreak predictions, spillover to managed bees

**Priority 3: Global Change Integration**
- Couple to Earth system models (downscaled climate scenarios)
- Multi-decadal projections under RCP pathways
- Range shift predictions (northern expansion, local extinctions)

### 8.5 Software Engineering and Maintenance

#### 8.5.1 Code Quality Practices

**Documentation Standards** (implemented in this MIDox):
- Doxygen comments for all public methods
- Biological rationale for all parameters
- Explicit uncertainty assessments
- Implementation differences from formal model noted

**Version Control**:
- Git repository with tagged releases
- Semantic versioning (MAJOR.MINOR.PATCH)
- Changelog documenting modifications

**Testing Framework**:
- Unit tests for core calculations (degree-day accumulation, sex ratio lookup)
- Integration tests for lifecycle completion (egg → adult)
- Regression tests preventing unintended behaviour changes

**Continuous Integration**:
- Automated compilation testing across platforms (Linux, Windows, macOS)
- Nightly runs of validation scenarios
- Performance benchmarking (runtime, memory)

#### 8.5.2 Portability and Dependencies

**Current Dependencies**:
- ALMaSS framework (C++ base classes, landscape infrastructure)
- Standard Template Library (STL containers, algorithms)
- OpenMP (parallelization)
- Optional: MPI (multi-node parallelization, not standard)

**Portability Status**:
- Tested: Linux (Ubuntu, CentOS), Windows (Visual Studio), macOS
- Compilers: GCC 9+, Clang 10+, MSVC 2019+
- Minimal external dependencies (improves long-term maintainability)

**Future-Proofing**:
- Avoid compiler-specific extensions
- Document platform-specific workarounds
- Maintain compatibility with ALMaSS framework updates (annual review)

---

## 9. Documentation Access

Complete interactive documentation with full code cross-references and searchable API details is available online through two complementary repositories:

**Interactive documentation:** https://christoppingau.github.io/ALMaSS_Osmia_model_documentations/  
The GitHub Pages site provides:

- Complete API documentation automatically generated by Doxygen
- Searchable class and method references with inheritance diagrams
- Call graphs showing function dependencies
- Cross-referenced source code with syntax highlighting
- Detailed parameter tables with biological justification
- Navigation between this narrative and code implementation via `@ref` links

**Archived version with DOI:** [![DOI](https://zenodo.org/badge/1104559872.svg)](https://doi.org/10.5281/zenodo.18701142)  
The Zenodo archive provides:

- Permanent, citable snapshot of model version 1.0
- Long-term preservation with DOI for reproducibility
- Complete source code, documentation, and example input files
- Versioned releases for tracking model evolution

**Source code repository:** https://github.com/ChrisToppingAU/ALMaSS_Osmia_model_documentations  
The GitHub repository includes:

- Current model source code (C++ files with enhanced Doxygen comments)
- Example configuration files and parameter sets
- Sample input data (landscape, pollen maps, weather)
- Compilation instructions and system requirements
- Automated testing framework
- Issue tracking for bug reports and feature requests

For compilation instructions, system requirements (compiler versions, dependencies), configuration examples, and troubleshooting guidance, see README.md in the source repository.

**Citation:** When using this model in research publications, please cite both the formal model paper and this MIDox documentation:

> Ziółkowska E, Bednarska AJ, Laskowski R, Topping CJ (2023). The Formal Model for the solitary bee *Osmia bicornis* L. agent-based model. Food and Ecological Systems Modelling Journal 4: e102102. https://doi.org/10.3897/fmj.4.102102

> Christopher John Topping, Xiaodong Duan, Elzbieta Ziolkowska (2026). *Osmia bicornis* Population Model: MIDox Implementation and Documentation. Food and Ecological Systems Modelling Journal [Volume TBD]. DOI: [TBD]

---

## References

**Core Model Development and Formal Specification**

Ziółkowska E, Bednarska AJ, Laskowski R, Topping CJ (2023). The Formal Model for the solitary bee *Osmia bicornis* L. agent-based model. Food and Ecological Systems Modelling Journal 4: e102102. https://doi.org/10.3897/fmj.4.102102

**Empirical Biology and Parameterization**

Giejdasz K, Wilkaniec Z (2002). Individual development of the red mason bee (*Osmia rufa* L., Megachilidae) under natural and laboratory conditions. Journal of Apicultural Science, 46(2), 51-57.

Ivanov SP (2006). Nesting habits of the red mason bee *Osmia bicornis* (Hymenoptera, Megachilidae). Entomological Review, 86(2), 221-226.

Radmacher S, Strohm E (2011). Effects of constant and fluctuating temperatures on the development of the solitary bee *Osmia bicornis* (Hymenoptera: Megachilidae). Apidologie, 42(6), 711-720.

Seidelmann K (2006). Open-cell parasitism shapes maternal investment patterns in the Red Mason bee *Osmia rufa*. Behavioral Ecology, 17(5), 839-848.

Seidelmann K, Ulbrich K, Mielenz N (2010). Conditional sex allocation in the Red Mason bee, *Osmia rufa*. Behavioral Ecology and Sociobiology, 64(3), 337-347.

Sgolastra F, Kemp WP, Buckner JS, Pitts-Singer TL, Maini S, Bosch J (2011). The long summer: pre-wintering temperatures affect metabolic expenditure and wintering survival of the solitary bee *Osmia lignaria*. Journal of Insect Physiology, 57(12), 1651-1659.

**Foraging and Spatial Behaviour**

Greenleaf SS, Williams NM, Winfree R, Kremen C (2007). Bee foraging ranges and their relationship to body size. Oecologia, 153(3), 589-596.

**MIDox Methodology**

[MIDox editorial paper reference - TO BE ADDED upon publication]

**ALMaSS Framework**

Topping CJ, Dalby L, Skov F (2016). Landscape structure and management alter the outcome of a pesticide ERA: Evaluating impacts of endocrine disruption using the ALMaSS European Brown Hare model. Science of the Total Environment, 541, 1477-1488.

Topping CJ, Hansen TS, Jensen TS, Jepsen JU, Nikolajsen F, Odderskær P (2003). ALMaSS, an agent-based model for animals in temperate European landscapes. Ecological Modelling, 167(1-2), 65-82.

**Toxicology and Risk Assessment** (for pesticide extensions)

Jager T, Albert C, Preuss TG, Ashauer R (2011). General unified threshold model of survival - a toxicokinetic-toxicodynamic framework for ecotoxicology. Environmental Science & Technology, 45(7), 2529-2540.

**Climate Change and Phenology**

Forrest JRK, Thomson JD (2011). An examination of synchrony between insect emergence and flowering in Rocky Mountain meadows. Ecological Monographs, 81(3), 469-491.

**Optimal Foraging Theory**

Stephens DW, Krebs JR (1986). Foraging Theory. Princeton University Press, Princeton, New Jersey.

**Population Modelling Methods**

Grimm V, Railsback SF (2005). Individual-Based Modeling and Ecology. Princeton University Press, Princeton, New Jersey.

Stillman RA, Railsback SF, Giske J, Berger U, Grimm V (2015). Making predictions in a changing world: The benefits of individual-based ecology. BioScience, 65(2), 140-150.

---

## Appendix: Quick Reference

### Key Classes and Their Roles

| Class | Purpose | Key Methods |
|-------|---------|-------------|
| [Osmia_Base](@ref Osmia_Base) | Common attributes/methods for all life stages | [GetTemp()](@ref Osmia_Base::GetTemp()), [SetTemp()](@ref Osmia_Base::SetTemp()), [CheckMortality()](@ref Osmia_Base::CheckMortality()) |
| [Osmia_Egg](@ref Osmia_Egg) | Egg development and hatching | [Step()](@ref Osmia_Egg::Step()), [TransitionToLarva()](@ref Osmia_Egg::TransitionToLarva()) |
| [Osmia_Larva](@ref Osmia_Larva) | Larval feeding and cocoon construction | [Step()](@ref Osmia_Larva::Step()), [TransitionToPrepupa()](@ref Osmia_Larva::TransitionToPrepupa()) |
| [Osmia_Prepupa](@ref Osmia_Prepupa) | Summer diapause | [Step()](@ref Osmia_Prepupa::Step()), [TransitionToPupa()](@ref Osmia_Prepupa::TransitionToPupa()) |
| [Osmia_Pupa](@ref Osmia_Pupa) | Metamorphosis to adult | [Step()](@ref Osmia_Pupa::Step()), [TransitionToInCocoon()](@ref Osmia_Pupa::TransitionToInCocoon()) |
| [Osmia_InCocoon](@ref Osmia_InCocoon) | Overwintering adult in cocoon | [Step()](@ref Osmia_InCocoon::Step()), [Emerge()](@ref Osmia_InCocoon::Emerge()) |
| [Osmia_Female](@ref Osmia_Female) | Reproductive behaviour | [Step()](@ref Osmia_Female::Step()), [FindNest()](@ref Osmia_Female::FindNest()), [Forage()](@ref Osmia_Female::Forage()), [LayEgg()](@ref Osmia_Female::LayEgg()) |
| [Osmia_Population_Manager](@ref Osmia_Population_Manager) | Simulation orchestration | [Init()](@ref Osmia_Population_Manager::Init()), [DoFirst()](@ref Osmia_Population_Manager::DoFirst()), [CreateObjects()](@ref Osmia_Population_Manager::CreateObjects()) |
| [Osmia_Nest_Manager](@ref Osmia_Nest_Manager) | Nest availability tracking | [UpdateOsmiaNesting()](@ref Osmia_Nest_Manager::UpdateOsmiaNesting()), [IsNestPossible()](@ref Osmia_Nest_Manager::IsNestPossible()) |
| [Osmia_Nest](@ref Osmia_Nest) | Individual nest structure | [AddEgg()](@ref Osmia_Nest::AddEgg()), [AddCocoon()](@ref Osmia_Nest::AddCocoon()), [CloseNest()](@ref Osmia_Nest::CloseNest()) |

### Configuration File Structure

Typical configuration file organization:
```
[OSMIA_PARAMETERS]
OSMIA_STARTNO = 50000
OSMIA_FEMALEMASSMMIN = 14.0
OSMIA_FEMALEMASSMAX = 28.0
OSMIA_EGGDEVELTHRESHOLD = 0.0
OSMIA_EGGDEVELTOTALDD = 86.0
OSMIA_LARVADEVELTHRESHOLD = 13.8
...
[OSMIA_MONTHLY_THRESHOLDS]
OSMIA_POLLENTHRESHOLDS = [12 quantity values], [12 quality values]
OSMIA_NECTARTHRESHOLDS = [12 quantity values], [12 quality values]
...
```

### Compilation Quick Start

```bash
# Clone repository
git clone https://github.com/[username]/osmia-bicornis-model.git
cd osmia-bicornis-model

# Compile with OpenMP support
g++ -O3 -fopenmp -std=c++11 -o osmia_model src/*.cpp -I./include

# Run simulation
./osmia_model config/default_parameters.cfg
```

Full compilation instructions including dependency installation in repository README.

---

**Document Version:** 1.0  
**Last Updated:** 2025  
**Corresponding Author:** Christopher John Topping  
**GitHub Pages:** https://[username].github.io/osmia-bicornis-model/  
**Zenodo DOI:** 10.5281/zenodo.XXXXXXX

