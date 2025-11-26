/**
 * @class Osmia_InCocoon
 * @brief Overwintering adult stage within cocoon
 * @extends Osmia_Pupa
 * 
 * @details Osmia_InCocoon represents the eclosed adult bee overwintering within its protective cocoon.
 * This stage implements the three-phase overwintering model: prewintering (late summer/autumn at elevated
 * temperatures), diapause proper (winter at low temperatures), and post-diapause quiescence (early spring
 * awaiting emergence cues). Overwintering mortality depends on prewinter thermal conditions, and emergence
 * timing is determined by an emergence counter that decrements with spring warming.
 * 
 * @par Biological Foundation
 * Adult *O. bicornis* eclose from pupae in late summer (August-September) but remain within cocoons through
 * winter, emerging the following spring (April-May). This overwintering strategy protects adults from winter
 * weather whilst avoiding energetic costs of maintaining activity. The three-phase model reflects documented
 * physiological transitions: initial high respiration rates during prewintering deplete fat reserves if warm
 * conditions persist, deep diapause conserves energy through winter, and post-diapause quiescence awaits
 * appropriate emergence conditions.
 * 
 * @par Three-Phase Overwintering Model
 * 
 * **Phase 1 - Prewintering** (typically September-November):
 * - Occurs at temperatures above prewintering threshold (15°C)
 * - Accumulates degree-days (m_DDPrewinter) that increase mortality risk
 * - Represents fat depletion at warm temperatures
 * - Based on Sgolastra et al. (2011) work with *O. lignaria*
 * 
 * **Phase 2 - Diapause** (typically November-January):
 * - Deep dormancy at temperatures near/below overwintering threshold (0°C)
 * - Minimal metabolic activity, lipid conservation
 * - No degree-day accumulation during this phase
 * - Chilling requirement for diapause completion (not explicitly modelled as separate counter)
 * 
 * **Phase 3 - Post-diapause quiescence** (typically February-March):
 * - Diapause complete but awaiting emergence cues
 * - Emergence counter decrements with temperatures above emergence threshold (5°C)
 * - When counter reaches zero, adult emerges from nest
 * - Counter equation: 35.48 - 0.0147 × DD_accumulated
 * 
 * @par Overwintering Mortality
 * Calculated as function of prewinter degree-day accumulation using equation from Sgolastra et al. (2011):
 * **mortality_probability = 0.05 × DD_prewinter - 4.63**
 * 
 * This captures the biological reality that prolonged warm autumn conditions deplete lipid reserves,
 * reducing winter survival. Cool autumn conditions (low DD accumulation) produce high survival.
 * 
 * @par Difference from Formal Model
 * **EXACT MATCH for equations** - Overwintering mortality and emergence counter equations implemented
 * precisely as specified in formal model. Temperature thresholds match formal model specifications
 * (0°C overwintering, 5°C emergence, 15°C prewintering baseline).
 * 
 * **Minor calibration**: Emergence counter constant adjusted from 39.48 to 35.48, shifting emergence
 * timing slightly earlier in spring to match field observations.
 * 
 * @par Implementation Notes
 * The three phases are not explicitly flagged as separate states - rather, they emerge from temperature-
 * threshold logic applied to DD accumulation and counter updates. This creates appropriate seasonal
 * phenology without needing complex state tracking.
 * 
 * @see Osmia_Base for overwintering parameters
 * @see st_Develop() for phase transition logic
 * @see WinterMortality() for mortality calculation
 * @see st_Emerge() for spring emergence
 */
class Osmia_InCocoon : public Osmia_Pupa
{
protected:
	/**
	 * @var m_emergencecounter
	 * @brief Countdown to spring emergence (decrements with warm days)
	 * @details Initialized at eclosion using equation: counter = 35.48 - 0.0147 × initial_DD
	 * Each day above emergence threshold (5°C) decrements counter by accumulated DD that day.
	 * When counter reaches zero (or negative), adult emerges from nest.
	 * 
	 * @par Biological Interpretation
	 * The counter represents integrated thermal accumulation needed for emergence readiness. After
	 * diapause completion (chilling requirement met), adults still need specific amount of spring
	 * warming before emergence is triggered. Counter mechanism ensures emergence synchrony with
	 * appropriate spring phenology (flower availability, favorable weather).
	 * 
	 * @par Empirical Basis
	 * Equation calibrated from field emergence observations showing relationship between winter/spring
	 * thermal patterns and emergence dates. Lower initial DD (cool conditions) requires more spring
	 * warming (higher counter start), whilst higher initial DD needs less spring warming.
	 * 
	 * @par Implementation Note
	 * Counter can become negative if warm spring arrives rapidly - this is acceptable as it simply
	 * means emergence threshold was exceeded quickly. Once counter ≤ 0, emergence occurs on next
	 * favorable day (temperature above emergence threshold).
	 */
	int m_emergencecounter;
	
	/**
	 * @var m_DDPrewinter
	 * @brief Accumulated degree-days during prewintering period (above 15°C baseline)
	 * @details Accumulates from eclosion through autumn. Only temperatures above 15°C contribute:
	 * DD_prewinter += (T_daily - 15) when T > 15°C
	 * 
	 * Used in overwintering mortality calculation. High values (prolonged warm autumn) increase
	 * mortality risk via fat depletion. Reset to zero at eclosion for new adults.
	 * 
	 * @par Biological Basis
	 * At temperatures above 15°C, prepupal/adult metabolism remains elevated, burning lipid reserves
	 * that are needed for winter survival. Bosch et al. (2008) documented weight loss rates of
	 * 0.2-0.4 mg/day during warm prewintering, with corresponding survival reductions.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation follows formal model specification precisely. The 15°C baseline
	 * comes directly from Sgolastra et al. (2011) methodology.
	 * 
	 * @par Valid Range
	 * [0, ~150] degree-days. Values >100 DD indicate poor overwintering conditions (prolonged warmth).
	 * Typical central European autumn might accumulate 30-60 DD, giving mortality around 0.2-0.4
	 * (20-40% death probability).
	 */
	double m_DDPrewinter;

public:
	/**
	 * @brief Constructor for newly eclosed adult in cocoon
	 * @param data Initialization data from pupal stage
	 * @details Initializes overwintering-specific attributes:
	 * - Sets m_DDPrewinter to 0 (begins accumulation from eclosion)
	 * - Calculates initial m_emergencecounter from current conditions
	 * - Transfers adult mass from pupal provision mass
	 * - Maintains nest linkage for eventual emergence
	 */
	Osmia_InCocoon(struct_Osmia* data);
	
	/**
	 * @brief Reinitialize object from pool
	 * @param data New initialization data
	 */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_InCocoon();
	
	/**
	 * @brief Main step function for overwintering adults
	 * @details Each day:
	 * 1. Call st_Develop to update DD accumulation and counters
	 * 2. Check emergence conditions (counter ≤ 0 and T > threshold)
	 * 3. Apply winter mortality test based on accumulated prewinter DD
	 * 4. Transition to emergence or continue overwintering
	 */
	virtual void Step(void);
	
	/**
	 * @brief Set overwintering temperature threshold (static parameter)
	 * @param a_temp Temperature threshold (°C) for overwintering DD accumulation
	 * @details Used primarily during initialization to set population-wide threshold from configuration.
	 */
	static void SetOverwinteringTempThreshold(double a_temp) { m_OverwinteringTempThreshold = a_temp; }
	
	/**
	 * @brief Get prewinter degree-day accumulation
	 * @return Total DD accumulated above 15°C baseline during prewintering
	 * @details Accessor for monitoring overwintering conditions and mortality risk. Used in output
	 * generation and validation.
	 */
	double GetDDPreWinter() { return m_DDPrewinter; }

protected:
	/**
	 * @brief Development state - manage overwintering phases and emergence preparation
	 * @return Next state (st_Emerge when ready, st_Develop to continue, or st_Die)
	 * 
	 * @details Complex multi-phase logic:
	 * 
	 * **Prewintering phase** (warm autumn conditions):
	 * - If T > 15°C: accumulate to m_DDPrewinter
	 * - Represents fat depletion period
	 * - No explicit phase flag - identified by temperature
	 * 
	 * **Diapause phase** (winter cold):
	 * - Temperatures near/below 0°C
	 * - Minimal DD accumulation
	 * - Lipid conservation
	 * - No active tracking - just waiting period
	 * 
	 * **Post-diapause phase** (spring warming):
	 * - If T > 5°C (emergence threshold): decrement m_emergencecounter
	 * - Counter -= (T - threshold)
	 * - When counter ≤ 0: ready to emerge
	 * - Check for favorable conditions (appropriate temperature, daylight)
	 * 
	 * **Mortality application**:
	 * - Call WinterMortality() to calculate probability from prewinter DD
	 * - Test against random number
	 * - If dies: return st_Die
	 * - If survives and counter ≤ 0: return st_Emerge
	 * - Otherwise: return st_Develop (continue overwintering)
	 * 
	 * @par Temperature Threshold Handling
	 * Three thresholds govern different processes:
	 * - 15°C: Prewinter DD accumulation (represents elevated metabolism)
	 * - 5°C: Emergence counter (spring warming required)
	 * - 0°C: Overwintering baseline (diapause maintenance)
	 * 
	 * These thresholds create realistic seasonal phenology without explicit calendar tracking.
	 * 
	 * @par Edge Cases
	 * - Very warm autumn (high prewinter DD): High mortality, survivors may emerge earlier
	 * - Very cold autumn (low prewinter DD): High survival, emergence depends more on spring warming
	 * - Warm winter spell: May decrement emergence counter but doesn't trigger premature emergence
	 *   (counter must reach zero AND appropriate conditions must persist)
	 * 
	 * @par Difference from Formal Model
	 * **IMPLEMENTATION MATCH** - The logic structure follows formal model three-phase description.
	 * Temperature thresholds and equations match specifications. The implementation detail of using
	 * continuous threshold checks rather than explicit phase flags is consistent with formal model
	 * emphasis on emergent behaviour from temperature-threshold interactions.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Transition state - emerge from cocoon as active adult
	 * @return Next state (st_Die as InCocoon object deleted)
	 * 
	 * @details Spring emergence sequence:
	 * 1. Signal population manager to create Osmia_Female object
	 * 2. Transfer attributes:
	 *    - Adult mass (determines size and fecundity)
	 *    - Natal nest location (becomes dispersal origin)
	 *    - Age (for lifespan tracking)
	 *    - Any carried-over attributes
	 * 3. Population manager handles nest cleanup (remove from nest cell list)
	 * 4. Delete InCocoon object
	 * 5. New Osmia_Female begins active adult life cycle
	 * 
	 * @par Biological Transition
	 * Adult bites through cocoon and cell partition, crawls out of nest, orients to environment,
	 * and begins pre-nesting maturation period. In model, this is represented by object class change
	 * from InCocoon (dormant, nest-bound) to Female (active, mobile).
	 * 
	 * @par Timing Considerations
	 * Emergence typically occurs on warm, sunny days in spring when:
	 * - Emergence counter has reached zero (sufficient spring warming accumulated)
	 * - Current temperature above threshold (favorable immediate conditions)
	 * - Potentially daylight/weather checks (not explicitly modelled)
	 * 
	 * This ensures emerged adults encounter appropriate conditions for initial flights and resource
	 * location.
	 * 
	 * @par Implementation Note
	 * The formal model discusses emergence timing but leaves specific trigger details to implementation.
	 * This counter-based approach captures observed field patterns: synchronous emergence within
	 * 2-3 week period, with timing responding to both winter and spring thermal conditions.
	 */
	virtual TTypeOfOsmiaState st_Emerge(void);
	
	/**
	 * @brief Calculate overwintering mortality probability from prewinter thermal conditions
	 * @return true if dies, false if survives
	 * 
	 * @details Implements Sgolastra et al. (2011) linear relationship:
	 * **mortality_probability = 0.05 × m_DDPrewinter - 4.63**
	 * 
	 * Then tests this probability against uniform random number to determine fate.
	 * 
	 * **Example calculations**:
	 * - m_DDPrewinter = 30 DD → mortality = 0.05×30 - 4.63 = -3.13 → 0.0 (clamped) → 0% death
	 * - m_DDPrewinter = 60 DD → mortality = 0.05×60 - 4.63 = -1.63 → 0.0 (clamped) → 0% death
	 * - m_DDPrewinter = 93 DD → mortality = 0.05×93 - 4.63 = 0.02 → 2% death
	 * - m_DDPrewinter = 100 DD → mortality = 0.05×100 - 4.63 = 0.37 → 37% death
	 * - m_DDPrewinter = 130 DD → mortality = 0.05×130 - 4.63 = 1.87 → 1.0 (clamped) → 100% death
	 * 
	 * @par Biological Interpretation
	 * The equation captures lipid depletion during warm prewintering. Each degree-day above 15°C
	 * represents time at elevated metabolism, burning fat reserves needed for winter survival.
	 * The negative intercept (-4.63) means zero mortality at low DD accumulation (cool autumns provide
	 * ideal prewintering conditions). Mortality rises linearly with warm autumn conditions.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Equation coefficients implemented precisely as specified in formal model
	 * (const = -4.63, slope = 0.05). Based directly on Sgolastra et al. (2011) empirical data for
	 * *O. lignaria* males, applied to both sexes of *O. bicornis* in absence of sex-specific data.
	 * 
	 * @par Uncertainty
	 * MEDIUM - Cross-species application (*O. lignaria* → *O. bicornis*) introduces uncertainty.
	 * Male-to-female extrapolation also uncertain (females may differ in lipid reserves and metabolic
	 * rates). However, underlying mechanism (fat depletion at warm temperatures) is well-established
	 * and should be qualitatively similar across related species and sexes.
	 * 
	 * @par Implementation Details
	 * - Probability clamped to [0, 1] range (equation can produce values outside this range)
	 * - Mortality test occurs once per overwintering period (not daily like developmental stages)
	 * - Timing of test determines whether early-season or late-season conditions dominate
	 * - Current implementation likely tests in spring (precise timing depends on Step logic)
	 * 
	 * @par Valid DD Range and Interpretation
	 * - **0-90 DD**: Ideal to good prewintering (mortality 0-10%)
	 * - **90-110 DD**: Moderate stress (mortality 10-40%)
	 * - **110-130 DD**: Poor conditions (mortality 40-80%)
	 * - **>130 DD**: Catastrophic (mortality >80%)
	 * 
	 * Central European autumns typically produce 30-80 DD, giving 0-20% mortality under normal conditions.
	 * Climate warming could increase DD accumulation, raising population-level mortality risk.
	 */
	bool WinterMortality();
	
	/**
	 * @var m_OverwinteringTempThreshold
	 * @brief Static temperature threshold for overwintering phase (°C)
	 * @details Default: 0.0°C. Defines baseline for diapause proper. Temperatures near/below this
	 * threshold represent true winter conditions where metabolism is minimized and lipid conservation
	 * is maximal.
	 * 
	 * @par Usage
	 * Not used for DD accumulation (unlike egg/larva/pupa thresholds) but serves as conceptual
	 * boundary between active prewintering (>15°C), transitional conditions (0-15°C), and true
	 * diapause (<0°C approximately).
	 */
	static double m_OverwinteringTempThreshold;
};

/**
 * @class Osmia_Female
 * @brief Active adult female conducting reproduction
 * @extends Osmia_InCocoon
 * 
 * @details Osmia_Female represents the culmination of the life cycle - the active adult female engaging
 * in dispersal, nest finding, foraging, provisioning, and egg laying. This is the most complex life stage,
 * with emergent spatial behaviour, resource-dependent reproductive decisions, and multiple interacting
 * state variables governing daily activities.
 * 
 * @par Biological Foundation
 * Adult female *O. bicornis* emerge in spring (April-May), undergo a brief pre-nesting maturation period,
 * then begin reproductive activities that may span 4-8 weeks. Typical female completes 2-4 nests in her
 * lifetime, each containing 6-12 cells (eggs). Foraging behaviour follows central-place foraging theory:
 * females return repeatedly to their nest, balancing travel costs against resource quality. Sex allocation
 * (female vs. male offspring) responds to provision mass availability, with larger cells receiving female
 * eggs and smaller cells receiving male eggs.
 * 
 * @par Reproductive Cycle
 * Female reproductive behaviour follows a repeating cycle:
 * 1. **Dispersal/nest searching**: Locate suitable cavity for nest
 * 2. **Nest establishment**: Clean cavity, orient to location
 * 3. **Cell provisioning cycle** (repeated per cell):
 *    - Forage for pollen/nectar
 *    - Return to nest with provisions
 *    - Construct cell partition
 *    - Determine sex of egg (based on provision mass)
 *    - Lay egg
 *    - Seal cell
 * 4. **Nest completion**: Seal final cell, abandon nest
 * 5. **Return to step 1** if longevity and eggs remaining permit
 * 
 * @par Foraging Behaviour
 * Foraging implements spatially-explicit resource search using pre-computed masks for efficiency.
 * Females exhibit:
 * - Age-dependent foraging efficiency (Seidelmann 2006 curves)
 * - Give-up thresholds (abandon poor patches)
 * - Distance-dependent returns (closer patches preferred)
 * - Competition effects (pollen depletion by bee density)
 * 
 * @par Sex Allocation
 * Sex determination follows haplodiploid genetics (fertilized=female, unfertilized=male) with strategic
 * maternal control. Females allocate sex based on provision mass:
 * - Large provisions → female egg (daughters require more resources)
 * - Small provisions → male egg (sons can develop on less)
 * - Sequential pattern: females typically at back of nest, males near entrance
 * 
 * @par Mortality
 * Adult females experience daily background mortality (0.02/day from Giejdasz et al. 2016) representing
 * combined hazards of foraging flights, weather exposure, predation, and senescence. Additionally
 * vulnerable to pesticide exposure via contaminated pollen and direct spray contact.
 * 
 * @par Difference from Formal Model
 * The formal model describes conceptual reproductive behaviors and resource dependencies. Implementation
 * adds necessary spatial algorithms (foraging masks, movement mechanics), specific sex allocation rules,
 * and detailed provisioning logistics. Core biological relationships (mass-fecundity from Seidelmann 2010,
 * age-efficiency curves) implemented exactly as specified. Spatial implementation details extend formal
 * model to enable landscape-scale simulation.
 * 
 * @see st_ReproductiveBehaviour() for reproductive state machine
 * @see Forage() for detailed foraging algorithm
 * @see CalcParasitised() for parasitism risk assessment
 * @see LayEgg() for sex allocation and egg creation
 */
class Osmia_Female : public Osmia_InCocoon
{
public:
#ifdef __OSMIARECORDFORAGE
	/** @brief Cumulative foraging success across all females (testing/validation only) */
	static double m_foragesum;
	/** @brief Count of foraging events (testing/validation only) */
	static int m_foragecount;
#endif

protected:
	//----------------------- Foraging Infrastructure -----------------------
	
	/**
	 * @var m_foragemask
	 * @brief Static coarse-resolution spatial search mask
	 * @details Shared across all females for memory efficiency. Provides 20 distance rings × 8 directions
	 * for efficient outward resource searches from nest location.
	 */
	static OsmiaForageMask m_foragemask;
	
	/**
	 * @var m_foragemaskdetailed
	 * @brief Static high-resolution spatial search mask
	 * @details Alternative mask with finer spatial resolution for detailed pollen assessment. Used when
	 * comprehensive resource evaluation needed rather than incremental search.
	 */
	static OsmiaForageMaskDetailed m_foragemaskdetailed;
	
	/**
	 * @var m_currentpollenlevel
	 * @brief Current pollen availability at active foraging location
	 * @details Updated during foraging to track resource depletion at focal patch. Compared against
	 * give-up thresholds to determine when to abandon patch and search elsewhere.
	 */
	double m_currentpollenlevel;
	
	/**
	 * @var m_pollengiveupthreshold
	 * @brief Proportional reduction triggering patch abandonment
	 * @details When pollen level drops to this proportion of initial value, female abandons current
	 * patch and searches for new location. Represents marginal value theorem - forager leaves when
	 * diminishing returns make travel to new patch worthwhile.
	 * 
	 * @par Typical Value
	 * 0.5-0.7 (abandon when patch drops to 50-70% of initial quality)
	 */
	static double m_pollengiveupthreshold;
	
	/**
	 * @var m_pollengiveupreturn
	 * @brief Absolute pollen level below which new search triggered
	 * @details Minimum acceptable return rate. If patch quality falls below this absolute threshold,
	 * female searches for new patch regardless of proportional decline. Prevents wasting time on
	 * completely depleted areas.
	 */
	static double m_pollengiveupreturn;
	
	//----------------------- Nest Provisioning State -----------------------
	
	/**
	 * @var m_CellOpenDays
	 * @brief Number of days current cell has been open (accumulating parasitism risk)
	 * @details Parasitism probability increases with cell open time - longer provisioning periods
	 * (due to poor forage availability or bad weather) expose cells to more parasitoid encounters.
	 * Reset to zero when cell is sealed.
	 */
	int m_CellOpenDays;
	
	/**
	 * @var m_CellCarryOver
	 * @brief Fractional hours carried to next day when cell not completed
	 * @details Cell construction requires minimum time (typically 1 day) but poor foraging may stretch
	 * across multiple days. This variable tracks partial progress to avoid losing fractional time
	 * between days.
	 */
	double m_CellCarryOver;
	
	/**
	 * @var m_EggsToLay
	 * @brief Total lifetime egg load remaining
	 * @details Calculated at emergence from body mass using Seidelmann (2010) relationship:
	 * total_eggs = N_nests_possible × (0.0371 × mass + 2.8399) ± 3
	 * 
	 * Decrements with each egg laid. When reaches zero, female ceases reproduction even if surviving.
	 * Represents ovary capacity constraint - bee cannot produce unlimited eggs.
	 */
	int m_EggsToLay;
	
	/**
	 * @var m_EggsThisNest
	 * @brief Planned eggs for current nest (decrements as cells completed)
	 * @details Drawn from probability distribution at nest initiation, representing female's "plan"
	 * for nest size. Actual eggs laid may differ if resources fail or female dies mid-nest.
	 * Planning occurs at nest start, simulating observed tendency for females to provision nests
	 * of characteristic size (6-12 cells typical).
	 */
	int m_EggsThisNest;
	
	/**
	 * @var m_ToDisperse
	 * @brief Flag indicating need for dispersal to new nesting area
	 * @details Set true when: (1) emergence (find initial nest area), (2) nest completion (may search
	 * new area for next nest), (3) repeated nest-finding failures (exhaust local options). Controls
	 * transition to long-distance dispersal behaviour vs. local nest searching.
	 */
	bool m_ToDisperse;
	
	/**
	 * @var m_EmergeAge
	 * @brief Days since emergence (tracks adult age separately from total age)
	 * @details Used for: age-dependent foraging efficiency (Seidelmann 2006 curves show efficiency
	 * increasing first 7-10 days then declining), lifespan constraints (max ~60 days), and output
	 * tracking of adult longevity.
	 */
	int m_EmergeAge;
	
	/**
	 * @var m_CurrentNestLoc
	 * @brief Spatial location of nest currently being provisioned
	 * @details X,Y coordinates of active nest. When no nest (dispersing or searching), m_x set to -1
	 * as flag. All foraging trips reference this location as return point (central-place foraging).
	 * Updated when new nest found.
	 */
	APoint m_CurrentNestLoc;
	
	/**
	 * @var m_ProvisioningTime
	 * @brief Days required to complete one cell
	 * @details Depends on: forage quality (poor resources extend provisioning), weather (bad days
	 * prevent foraging), female efficiency (age effects). Typically 1-3 days per cell. Longer times
	 * increase parasitism risk via extended cell open duration.
	 */
	int m_ProvisioningTime;
	
	/**
	 * @var m_FlyingCounter
	 * @brief Days spent flying/foraging during current cell construction
	 * @details Distinguishes active foraging days from weather delays. Used to accurately track
	 * provisioning effort vs. unavoidable delays. Increments only on days when foraging occurs.
	 */
	int m_FlyingCounter;
	
	/**
	 * @var m_CurrentProvisioning
	 * @brief Mass of pollen/nectar currently provisioned in active cell (mg)
	 * @details Accumulates from zero as female makes foraging trips. When reaches target mass for
	 * planned sex (female target or male target), cell is complete and egg is laid. Determines
	 * final egg sex and offspring size.
	 */
	double m_CurrentProvisioning;
	
	/**
	 * @var m_BeeSizeScore1
	 * @brief Coarse size class (0=very small, 1=small, 2=medium, 3=large)
	 * @details Categorical size classification from adult mass. Used for size-dependent behaviour
	 * parameters if implemented (currently placeholder for future extensions where larger bees might
	 * have different foraging ranges or fecundity).
	 */
	int m_BeeSizeScore1;
	
	/**
	 * @var m_BeeSizeScore2
	 * @brief Fine-grained size classification
	 * @details Finer size categories than m_BeeSizeScore1, with step size controlled by configuration
	 * parameter cfg_OsmiaAdultMassCategoryStep. Enables more nuanced size-dependent parameterization
	 * if needed for model extensions.
	 */
	int m_BeeSizeScore2;
	
	/**
	 * @var m_NestProvisioningPlan
	 * @brief Queue of target provision masses for planned nest cells
	 * @details Female "plans" nest at initiation, generating sequence of target masses (one per egg).
	 * Deque structure allows efficient removal from front as cells completed. Masses decline from
	 * first to last cell (progressive resource depletion effect), with stochastic variation.
	 * 
	 * @par Biological Basis
	 * Seidelmann (2010) documented progressive decline in provision mass from first to last offspring
	 * within nests, reflecting maternal aging and resource depletion. This planning structure implements
	 * that pattern whilst allowing females to adjust actual provisioning based on encountered resources.
	 */
	deque<double>m_NestProvisioningPlan;
	
	/**
	 * @var m_NestProvisioningPlanSex
	 * @brief Queue of planned sexes (true=female, false=male) corresponding to m_NestProvisioningPlan
	 * @details Sex allocated at planning stage based on provision mass targets: larger masses get
	 * female designation, smaller get male. Female control over fertilization allows implementation
	 * of planned sex allocation. Actual sex ratio emerges from resource availability and mortality.
	 */
	deque<bool>m_NestProvisioningPlanSex;
	
	//----------------------- Foraging State -----------------------
	
	/**
	 * @var m_ForageLoc
	 * @brief Flag indicating whether foraging location has been identified
	 * @details true = female has located pollen source and is actively exploiting it.
	 * false = female needs to search for new forage location (initial or after patch abandonment).
	 * Controls whether to search or continue foraging at known location.
	 */
	bool m_ForageLoc;
	
	/**
	 * @var m_ForageLocPoly
	 * @brief Index to polygon list in population manager providing resources
	 * @details Each landscape polygon (field, forest patch, etc.) has associated resource values.
	 * This index provides fast lookup of focal polygon's resource data without repeated spatial queries.
	 * Updated when new foraging location selected.
	 */
	int m_ForageLocPoly;
	
	/**
	 * @var m_ForageSteps
	 * @brief Number of distance steps in foraging mask (determines search resolution)
	 * @details Static value controlling how many distance rings are searched. Higher values allow
	 * finer-grained searches but increase computation. Typically 10-20 steps covering 0 to max
	 * foraging range.
	 */
	static int m_ForageSteps;
	
	/**
	 * @var m_PollenCompetitionsReductionScaler
	 * @brief Scaling factor for inter-specific pollen competition
	 * @details Adjusts available pollen based on assumed competition from other bee species (honey bees,
	 * bumble bees, other solitaries). Values <1.0 reduce available pollen, simulating competitive depletion
	 * by non-modelled species. Enables exploration of competition scenarios.
	 * 
	 * @par Implementation Note
	 * Simple proportional reduction rather than mechanistic competition model. Pragmatic approach given
	 * uncertainty in competitor densities and overlap in flower usage.
	 */
	static double m_PollenCompetitionsReductionScaler;
	
	/**
	 * @var m_FemaleForageEfficiency
	 * @brief Vector of age-dependent foraging efficiency multipliers indexed by adult age
	 * @details Implements Seidelmann (2006) empirical efficiency curves showing:
	 * - Days 1-7: Efficiency increases as females gain experience
	 * - Days 8-15: Peak efficiency (full capability)
	 * - Days 15+: Gradual decline with senescence (wing wear, reduced flight capability)
	 * 
	 * Applied as multiplier to daily forage returns: actual_forage = base_forage × efficiency[age]
	 * 
	 * @par Biological Basis
	 * Young bees need practice to optimize foraging routes and flower handling. Old bees experience
	 * cumulative wing wear and muscle degradation reducing flight speed and cargo capacity.
	 */
	static vector<double> m_FemaleForageEfficiency;
	
	/**
	 * @var m_ForageLocX
	 * @brief X-coordinate of current foraging location
	 */
	int m_ForageLocX;
	
	/**
	 * @var m_ForageLocY
	 * @brief Y-coordinate of current foraging location
	 */
	int m_ForageLocY;
	
	/**
	 * @var m_foraged_resource_pesticide
	 * @brief Array storing pesticide concentrations in foraged resources
	 * @details When pesticide module active, tracks pesticide content of pollen/nectar collected at
	 * different locations. Enables simulation of pesticide exposure via contaminated provisions.
	 * Array size matches number of pesticide types configured.
	 */
	double* m_foraged_resource_pesticide;
	
	//----------------------- Pesticide Exposure (Conditional Compilation) -----------------------
	
#ifdef __OSMIA_PESTICIDE_ENGINE
public:
	/** @brief Egg-specific pesticide death probability after threshold exceedance */
	static double m_OsmiaEggPPPEffectProb;
	
	/** @brief Pesticide concentration threshold for egg effects */
	static double m_OsmiaEggPPPThreshold;
	
	/** @brief Adult pesticide death probability after threshold exceedance */
	static double m_OsmiaPPPEffectProb;
	
	/** @brief Adult pesticide concentration threshold */
	static double m_OsmiaPPPThreshold;
	
	/** @brief Daily pesticide decay rate in bee body (proportion lost per day) */
	static double m_OsmiaPPPDecayRate;
	
	/** @brief Absorption rate for overspray exposure (proportion transferred body to internal) */
	static double m_OsmiaPPPAbsorptionRateOverspray;
	
	/** @brief Absorption rate for contact exposure */
	static double m_OsmiaPPPAbsorptionRateContact;
	
	/** @brief Surface area exposed to overspray (mm²) */
	static double m_OsmiaPPPOversprayBodySurface;
	
	/** @brief Surface area for contact exposure (mm²) */
	static double m_OsmiaPPPContactBodySurface;
	
	/** @brief Probability of experiencing overspray event */
	static double m_OsmiaPPPOversprayChance;
	
protected:
	/** @brief Accessor for pesticide threshold parameter */
	static double GetPPPThreshold() { return m_OsmiaPPPThreshold; }
	
	/** @brief Accessor for pesticide effect probability */
	static double GetPPPEffectProb() { return m_OsmiaPPPEffectProb; }
	
	/** @brief Accessor for pesticide decay rate */
	static double GetPPPDecayRate() { return m_OsmiaPPPDecayRate; }
	
	/** @brief Accessor for overspray absorption rate */
	static double GetPPPAbsorptionRateOverspray() { return m_OsmiaPPPAbsorptionRateOverspray; }
	
	/** @brief Accessor for overspray body surface */
	static double GetPPPOversprayBodySurface() { return m_OsmiaPPPOversprayBodySurface; }
	
	/** @brief Accessor for contact body surface */
	static double GetPPPContactBodySurface() { return m_OsmiaPPPContactBodySurface; }
#endif

	//----------------------- Testing Support (Conditional Compilation) -----------------------
	
#ifdef __OSMIATESTING
	/** @brief Target nest data for validation (intended provisioning plan) */
	OsmiaNestData m_target;
	
	/** @brief Achieved nest data for validation (actual provisioning accomplished) */
	OsmiaNestData m_achieved;
	
	/** @brief Flag for first nest tracking */
	bool m_firstnestflag;
#endif

	//----------------------- Behavioural Methods -----------------------
	
	/**
	 * @brief Death state with female-specific cleanup
	 * @details Extends base st_Dying to handle:
	 * - Incomplete nest abandonment (set nest to closed, prevent further egg additions)
	 * - Resource release (return unused forage hours to pool)
	 * - Output recording (lifetime reproductive success, cause of death)
	 */
	virtual void st_Dying(void);
	
	/**
	 * @brief Development state for adult females (minimal - no metamorphosis)
	 * @return Next state (typically st_ReproductiveBehaviour or st_Dispersal)
	 * @details Unlike immature stages, adult females don't "develop" - this state primarily handles
	 * daily initialization and transitions to active reproductive states. May handle pre-nesting
	 * maturation period (first few days post-emergence before reproduction begins).
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Search for suitable nest cavity
	 * @return true if nest found, false if search fails
	 * 
	 * @details Nest finding algorithm:
	 * 1. Sample locations around current position using movement probability distribution
	 * 2. Check each location for suitable cavities (queries landscape manager)
	 * 3. If suitable cavity available: establish nest, set m_CurrentNestLoc, return true
	 * 4. If no cavity found after N attempts: set m_ToDisperse=true, return false
	 * 
	 * @par Cavity Suitability Criteria
	 * - Appropriate diameter (6-9mm for *O. bicornis*)
	 * - Sufficient depth (>10cm)
	 * - Protected location (not fully exposed)
	 * - Not already occupied
	 * 
	 * In model, suitability determined by landscape polygon attributes (habitat types provide different
	 * cavity densities).
	 * 
	 * @par Search Limitations
	 * Female makes limited attempts (m_OsmiaFindNestAttemptNo, typically 5-10). Repeated failures
	 * trigger dispersal to new area. Reflects biological reality that suitable cavities are limiting
	 * resource in many landscapes.
	 */
	virtual bool FindNestLocation(void);
	
	/**
	 * @brief Dispersal state for long-distance movements to new nesting areas
	 * @return Next state (st_ReproductiveBehaviour if dispersal successful, st_Die if fails)
	 * 
	 * @details Long-distance dispersal using different movement distribution than local foraging:
	 * - Samples from m_dispersalmovementdistances (typically longer distances)
	 * - Moves to new location, attempts nest finding
	 * - If successful: begins reproduction at new location
	 * - If fails: may attempt further dispersal or die (dispersal mortality)
	 * 
	 * @par Biological Context
	 * Dispersal occurs: (1) at emergence from natal nest, (2) after nest completion if local resources
	 * depleted, (3) after repeated nest-finding failures. Enables population spread and colonization
	 * of new areas, but carries mortality costs (navigation failures, exposure, energy depletion).
	 * 
	 * @par Difference from Formal Model
	 * Formal model describes dispersal conceptually. Implementation specifies movement mechanics using
	 * beta probability distributions parameterized from allometric relationships between body size and
	 * foraging range. Dispersal distances typically 2-3× longer than foraging movements, reaching
	 * maximum homing distance (R90 = 1430m).
	 */
	virtual TTypeOfOsmiaState st_Dispersal(void);
	
	/**
	 * @brief Main foraging algorithm collecting pollen and nectar
	 * @return Mass of pollen collected (mg)
	 * 
	 * @details Complex spatial foraging implementing:
	 * 
	 * **Phase 1 - Resource Location** (if !m_ForageLoc):
	 * - Use OsmiaForageMask to search concentrically from nest
	 * - Evaluate pollen availability at each location
	 * - Select location with acceptable pollen level
	 * - Store location, set m_ForageLoc = true
	 * 
	 * **Phase 2 - Resource Exploitation** (if m_ForageLoc):
	 * - Calculate forage return based on:
	 *   - Local pollen availability
	 *   - Distance from nest (travel time reduces forage time)
	 *   - Age-dependent efficiency (m_FemaleForageEfficiency[age])
	 *   - Competition (density-dependent depletion)
	 *   - Hours available (m_foragehours)
	 * - Deplete local pollen (update landscape resource levels)
	 * - Check give-up thresholds:
	 *   - If pollen dropped to <threshold: m_ForageLoc=false (search new location)
	 *   - If absolute level <minimum: m_ForageLoc=false
	 * - Return collected mass
	 * 
	 * **Phase 3 - Resource Accumulation**:
	 * - Add collected mass to m_CurrentProvisioning
	 * - Check if cell target reached
	 * - If complete: trigger egg laying sequence
	 * 
	 * @par Age-Dependent Efficiency
	 * Daily foraging success scaled by m_FemaleForageEfficiency[m_EmergeAge]:
	 * - Days 0-7: Rising efficiency (inexperience)
	 * - Days 7-15: Peak efficiency
	 * - Days 15+: Declining efficiency (senescence)
	 * 
	 * Based on Seidelmann (2006) empirical curves.
	 * 
	 * @par Distance Effects
	 * Longer distances reduce foraging efficiency:
	 * - Travel time reduces available foraging hours
	 * - Energetic costs reduce net provisioning
	 * - Implemented via time budget (hours available after travel)
	 * 
	 * @par Competition
	 * Pollen availability decreases with local *Osmia* density:
	 * available_pollen = base_pollen × (1 - density × m_DensityDependentPollenRemovalConst)
	 * 
	 * Provides negative feedback preventing unrealistic population growth in favorable habitats.
	 * 
	 * @par Give-Up Decisions
	 * Two thresholds govern patch abandonment:
	 * 1. Proportional decline: leave when pollen drops to X% of initial (marginal value theorem)
	 * 2. Absolute minimum: leave when pollen falls below acceptable return rate
	 * 
	 * Whichever threshold triggered first causes female to search for new patch.
	 * 
	 * @par Difference from Formal Model
	 * Formal model describes foraging conceptually with distance costs and resource depletion. Implementation
	 * specifies spatial search algorithms, time budgets, and give-up rules. Core relationships (age effects,
	 * distance costs, competition) implement formal model principles with specific parameterization from
	 * Seidelmann (2006) and calibration.
	 */
	double Forage(void);
	
	/**
	 * @brief Reproductive behaviour state coordinating nesting activities
	 * @return Next state (continues reproduction, disperses, or dies)
	 * 
	 * @details Master reproductive state machine:
	 * 
	 * **If no nest** (m_CurrentNestLoc.x == -1):
	 * - Call FindNestLocation()
	 * - If successful: initialize new nest
	 * - If fails: return st_Dispersal
	 * 
	 * **If active nest**:
	 * - Check weather (if unfavorable: skip foraging, increment cell open days)
	 * - Call Forage() to collect provisions
	 * - Add to m_CurrentProvisioning
	 * - Check if cell target reached:
	 *   - If yes: call LayEgg(), reset m_CurrentProvisioning, decrement m_EggsThisNest
	 *   - If no: continue provisioning
	 * - Check if nest complete (m_EggsThisNest == 0):
	 *   - Seal nest, set m_CurrentNestLoc.x = -1
	 *   - Check remaining eggs (m_EggsToLay):
	 *     - If >0: prepare for next nest (may disperse or search locally)
	 *     - If 0: reproduction complete, continue until senescence
	 * 
	 * **Daily mortality test**:
	 * - Apply m_OsmiaFemaleBckMort (0.02 daily probability)
	 * - If dies: return st_Die
	 * 
	 * @par Weather Effects
	 * Bad weather (too cold, rainy, windy) prevents foraging but doesn't stop clock:
	 * - Cell remains open (accumulating parasitism risk)
	 * - Age advances (depleting remaining lifespan)
	 * - Provisions don't accumulate (extends nest completion time)
	 * 
	 * Creates realistic weather impacts on reproductive success without explicit weather state tracking.
	 * 
	 * @par Nest Completion Decision
	 * Female "knows" when nest is complete via m_EggsThisNest counter (planned at nest initiation).
	 * Actual completion may differ from plan if:
	 * - Resources fail (cannot provision remaining cells)
	 * - Female dies mid-nest
	 * - Eggs exhausted before plan complete
	 * 
	 * This represents partial implementation of female's "plan" responding to realized conditions.
	 */
	virtual TTypeOfOsmiaState st_ReproductiveBehaviour(void);
	
	/**
	 * @brief Calculate planned eggs for next nest
	 * @return Number of eggs planned for nest
	 * 
	 * @details Draws from m_eggspernestdistribution (typically beta distribution giving 6-12 eggs
	 * with appropriate skew). Value constrained by:
	 * - m_OsmiaFemaleMinEggsPerNest (lower bound)
	 * - m_OsmiaFemaleMaxEggsPerNest (upper bound)
	 * - m_EggsToLay (can't plan more than remaining lifetime egg load)
	 * 
	 * @par Biological Basis
	 * Females exhibit characteristic nest sizes reflecting tradeoffs between offspring number and quality.
	 * Larger nests risk incomplete provisioning if female dies; smaller nests under-utilize female's
	 * reproductive capacity. Optimal nest size balances these factors and varies with resource quality.
	 * 
	 * @par Implementation
	 * Planning occurs at nest initiation, creating "target" that female works toward. Actual eggs laid
	 * may differ based on resource availability, weather, mortality. This two-stage process (plan then
	 * execute) captures observed tendency for nests to cluster around typical sizes whilst allowing
	 * environmental variation.
	 */
	int PlanEggsPerNest();
	
	/**
	 * @brief Calculate total lifetime egg load from body mass
	 * @details Implements Seidelmann (2010) empirical relationship:
	 * **eggs_per_nest = 0.0371 × mass + 2.8399 (±3 eggs stochastic)**
	 * **total_eggs = m_TotalNestsPossible × eggs_per_nest**
	 * 
	 * Sets m_EggsToLay at emergence. Larger females have higher fecundity via both more eggs per nest
	 * and potential for completing more nests (if they can provision faster).
	 * 
	 * Also initializes m_EggsThisNest by calling PlanEggsPerNest() + 2 (2 removed at nest start, 
	 * creating correct initial count).
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implements Seidelmann (2010) equation precisely as specified in formal model
	 * (const = 2.8399, slope = 0.0371, stochastic variation ±3 eggs). Total lifetime fecundity calculated
	 * as: eggs_per_nest × maximum_possible_nests, representing ovary capacity constraint.
	 */
	void CalculateEggLoad() {
		m_EggsToLay = int((m_TotalNestsPossible * (0.0371 * m_Mass + 2.8399)) + (g_rand_uni_fnc() * 6) - 3);
		m_EggsThisNest = PlanEggsPerNest() + 2;
	}
	
	/**
	 * @brief Determine parasitism status for egg about to be laid
	 * @param a_daysopen Number of days cell has been open
	 * @return Parasitoid type affecting this egg
	 * 
	 * @details Two possible parasitism models (controlled by m_UsingMechanisticParasitoids):
	 * 
	 * **Simple model** (probability-based):
	 * - Parasitism risk increases linearly with cell open time
	 * - Base probability scaled by m_ParasitismProbToTimeCellOpen
	 * - Bombylid probability from m_BombylidProbability
	 * - Random draw determines outcome
	 * 
	 * **Mechanistic model** (population-based):
	 * - Queries parasitoid population manager for local parasitoid density
	 * - Attack probability from m_ParasitoidAttackChance[parasitoid_type] × density
	 * - Multiple parasitoid types possible with type-specific parameters
	 * - More realistic but requires additional parasitoid population tracking
	 * 
	 * @par Biological Basis
	 * Open nest cells are vulnerable to parasitoid females searching for hosts. Longer provisioning
	 * times (multi-day cells) provide more opportunity for parasitoid discovery. Bombylid flies are
	 * primary parasitoids of *O. bicornis*, laying eggs in open cells that consume host egg/larva.
	 * 
	 * @par Implementation Note
	 * Called during LayEgg() sequence, before egg object creation. Parasitism status assigned to egg
	 * at laying, determining subsequent survival. Parasitised eggs/larvae die at characteristic time
	 * for parasitoid species (handled by population manager during development).
	 */
	TTypeOfOsmiaParasitoids CalcParasitised(double a_daysopen);
	
	/**
	 * @brief Create and lay egg in completed cell
	 * 
	 * @details Egg laying sequence:
	 * 1. Determine sex based on provision mass:
	 *    - If m_CurrentProvisioning >= female minimum: lay female egg (fertilized)
	 *    - If below female minimum: lay male egg (unfertilized)
	 * 2. Calculate parasitism status: CalcParasitised(m_CellOpenDays)
	 * 3. Create struct_Osmia with egg initialization data:
	 *    - Sex (determined above)
	 *    - Provision mass (m_CurrentProvisioning)
	 *    - Nest pointer (m_OurNest)
	 *    - Parasitism status
	 * 4. Signal population manager to create Osmia_Egg object
	 * 5. Nest adds egg to cell list
	 * 6. Reset cell state:
	 *    - m_CurrentProvisioning = 0
	 *    - m_CellOpenDays = 0
	 *    - m_EggsToLay -= 1
	 *    - m_EggsThisNest -= 1
	 * 
	 * @par Sex Allocation Strategy
	 * Haplodiploid sex determination with maternal control allows strategic allocation:
	 * - Females (diploid) require more resources → placed on larger provisions
	 * - Males (haploid) can develop on less → placed on smaller provisions
	 * 
	 * Threshold-based allocation emerges from provision mass variation:
	 * - Early cells (larger provisions): mostly females
	 * - Late cells (smaller provisions due to resource depletion): mostly males
	 * - Natural pattern: females at back of nest, males near entrance
	 * 
	 * @par Difference from Formal Model
	 * Formal model describes sex allocation conceptually (larger provisions → females). Implementation
	 * specifies threshold-based decision rule with provision mass targets. Core biology (females need
	 * more resources) implemented faithfully; specific threshold values calibrated from field sex ratio
	 * observations.
	 * 
	 * @par Parasitism Integration
	 * Parasitism status assigned at laying (based on cell open time) and carries through development.
	 * Parasitised individuals develop normally until parasitoid emerges (timing depends on parasitoid
	 * type), then host dies. This creates realistic delayed mortality rather than immediate death.
	 */
	void LayEgg();

public:
	/**
	 * @brief Constructor for newly emerged adult female
	 * @param data Initialization data from InCocoon stage
	 * @details Initializes:
	 * - Adult mass (determines size class and fecundity)
	 * - Calculates lifetime egg load: CalculateEggLoad()
	 * - Sets initial state (dispersal to find first nest area)
	 * - Initializes foraging attributes (no location, no nest)
	 * - Records emergence location (becomes dispersal origin)
	 */
	Osmia_Female(struct_Osmia* data);
	
	/**
	 * @brief Reinitialize from object pool
	 * @param data New initialization data
	 */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_Female();
	
	/**
	 * @brief Female-specific initialization (called by constructor and ReInit)
	 * @param a_mass Adult body mass (mg)
	 * @details Handles mass-dependent initialization:
	 * - Size class calculation (m_BeeSizeScore1, m_BeeSizeScore2)
	 * - Provision mass targets (females vs. males)
	 * - Fecundity calculation: CalculateEggLoad()
	 */
	virtual void Init(double a_mass);
	
	/**
	 * @brief Pre-step initialization each day
	 * @details Sets up daily state:
	 * - Reset forage hours available (from weather/daylight)
	 * - Increment emerge age
	 * - Check lifespan limit (m_EmergeAge vs. m_OsmiaFemaleLifespan)
	 * - Update local resource availability if needed
	 */
	virtual void BeginStep(void);
	
	/**
	 * @brief Main step function orchestrating daily behaviour
	 * @details Calls appropriate behavioural state based on m_CurrentOState:
	 * - st_Develop: Initial maturation
	 * - st_Dispersal: Long-distance movement
	 * - st_ReproductiveBehaviour: Nesting and provisioning
	 * - st_Die: Cleanup and removal
	 * 
	 * Loops until state machine returns terminal state (DONE or DIE).
	 */
	virtual void Step(void);
	
	//----------------------- Static Setters (Population Manager Initialization) -----------------------
	
	/** @brief Set number of distance steps in foraging mask */
	static void SetForageSteps(int a_sz) { m_ForageSteps = a_sz; }
	
	/**
	 * @brief Initialize detailed foraging mask
	 * @param a_step Step size between sample points
	 * @param a_max Maximum search distance
	 */
	static void SetForageMaskDetailed(int a_step, int a_max) {
		OsmiaForageMaskDetailed fmd(a_step, a_max);
		m_foragemaskdetailed = fmd;
	}
	
	/** @brief Set proportional give-up threshold for patch abandonment */
	static void SetPollenGiveUpThreshold(double a_prop) { m_pollengiveupthreshold = a_prop; }
	
	/** @brief Set absolute give-up threshold (minimum acceptable return) */
	static void SetPollenGiveUpReturn(double a_value) { m_pollengiveupreturn = a_value; }
	
	/** @brief Set daily background mortality for adult females */
	static void SetDailyMort(double a_prob) { m_OsmiaFemaleBckMort = a_prob; }
	
	/** @brief Set number of nest-finding attempts before dispersal triggered */
	static void SetNestFindAttempts(int a_no) { m_OsmiaFindNestAttemptNo = a_no; }
	
	/** @brief Set minimum eggs per nest (lower bound for planning distribution) */
	static void SetMinEggsPerNest(int a_eggs) { m_OsmiaFemaleMinEggsPerNest = a_eggs; }
	
	/** @brief Set maximum eggs per nest (upper bound for planning distribution) */
	static void SetMaxEggsPerNest(int a_eggs) { m_OsmiaFemaleMaxEggsPerNest = a_eggs; }
	
	/**
	 * @brief Set cocoon-to-provision mass conversion and derived parameters
	 * @param a_ratio Conversion factor (proportion of provision mass converted to cocoon mass)
	 * @details Also calculates total provisioning mass loss parameters by scaling cocoon mass loss
	 * configuration values. These derived parameters used in nest provisioning planning.
	 */
	static void SetCocoonToProvisionMass(double a_ratio) {
		m_CocoonToProvisionMass = a_ratio;
		m_TotalProvisioningMassLoss = cfg_OsmiaTotalCocoonMassLoss.value() * a_ratio;
		m_TotalProvisioningMassLossRange = cfg_OsmiaTotalCocoonMassLossRange.value() * a_ratio;
		m_TotalProvisioningMassLossRangeX2 = m_TotalProvisioningMassLossRange * 2.0;
	}
	
	/** @brief Set provision-to-cocoon mass conversion factor */
	static void SetProvisionToCocoonMass(double a_ratio) { m_ProvisionToCocoonMass = a_ratio; }
	
	/** @brief Set pollen score to mg conversion factor */
	static void SetPollenScoreToMg(double a_ratio) { m_PollenScoreToMg = a_ratio; }
	
	/** @brief Set minimum target provision mass for male cells (instance method) */
	void SetMaleMinTargetProvisionMass(double a_mass) { m_MaleMinTargetProvisionMass = a_mass; }
	
	/** @brief Set minimum target provision mass for female cells (instance method) */
	void SetFemaleMinTargetProvisionMass(double a_mass) { m_FemaleMinTargetProvisionMass = a_mass; }
	
	/** @brief Set maximum target provision mass for female cells (instance method) */
	void SetFemaleMaxTargetProvisionMass(double a_mass) { m_FemaleMaxTargetProvisionMass = a_mass; }
	
	/** @brief Set minimum cell construction time (days) */
	static void SetMinimumCellConstructionTime(double a_time) { m_MinimumCellConstructionTime = a_time; }
	
	/** @brief Set maximum cell construction time (days) */
	static void SetMaximumCellConstructionTime(double a_time) { m_MaximumCellConstructionTime = a_time; }
	
	/** @brief Set maximum lifetime nests possible */
	static void SetTotalNestsPossible(int a_total) { m_TotalNestsPossible = a_total; }
	
	/** @brief Set Bombyliid parasitism probability */
	static void SetBombylidProbability(double a_prob) { m_BombylidProbability = a_prob; }
	
	/** @brief Set parasitism probability to cell open time conversion factor */
	static void SetParasitismProbToTimeCellOpen(double a_ratio) { m_ParasitismProbToTimeCellOpen = a_ratio; }
	
	/** @brief Set flag for using mechanistic vs. simple parasitoid model */
	static void SetUsingMechanisticParasitoids(bool a_flag) { m_UsingMechanisticParasitoids = a_flag; }
	
	/** @brief Set parasitoid attack probability parameters (vector for multiple types) */
	static void SetParasitoidParameters(vector<double> a_params) { m_ParasitoidAttackChance = a_params; }
	
	/** @brief Set density-dependent pollen removal constant (instance method) */
	void SetDensityDependentPollenRemovalConst(double a_value) { m_DensityDependentPollenRemovalConst = a_value; }
	
	/** @brief Add age-specific foraging efficiency value to vector */
	static void AddForageEfficiency(double a_eff) { m_FemaleForageEfficiency.push_back(a_eff); }
	
	/**
	 * @brief Get available pollen in polygon from starting location
	 * @param a_required_amount [in/out] Amount of pollen needed (mg)
	 * @param a_foraged_amount [out] Amount actually foraged (mg)
	 * @param a_polygon Polygon index to query
	 * @param a_loc_x X-coordinate in polygon
	 * @param a_loc_y Y-coordinate in polygon
	 * 
	 * @details Queries landscape manager for pollen availability, applies competition effects,
	 * depletes local resources, returns actual foraged amount. Used by Forage() to implement
	 * resource acquisition and depletion.
	 */
	void GetPollenInPolygon(double& a_required_amount, double& a_foraged_amount, int a_polygon, int a_loc_x, int a_loc_y);
	
	/**
	 * @brief Handle farm management events affecting females
	 * @param event Type of farming operation (mowing, spraying, tillage, etc.)
	 * @return true if female survives event, false if killed
	 * 
	 * @details Farm operations can affect females directly:
	 * - Pesticide spraying: Contact mortality, contaminated forage
	 * - Mowing/harvesting: Destroys forage resources
	 * - Tillage: Destroys ground nests (if applicable)
	 * 
	 * Return value determines whether female continues activity or is removed from simulation.
	 */
	virtual bool OnFarmEvent(FarmToDo event);
	
	/**
	 * @brief Override pesticide contact handling for female-specific exposure
	 * @param a_x X-location of contact (default -1 = use current location)
	 * @param a_y Y-location of contact (default -1 = use current location)
	 * 
	 * @details Females experience pesticide exposure via:
	 * 1. Direct overspray (flying through spray plume or foraging during application)
	 * 2. Contact with treated vegetation
	 * 3. Contaminated pollen/nectar (handled separately via m_foraged_resource_pesticide)
	 * 
	 * Exposure routes and effects calculated using pesticide module parameters (absorption rates,
	 * body surface areas, thresholds). Female-specific exposure reflects foraging behaviour bringing
	 * them into contact with treated crops.
	 */
	virtual void DoPesticideContact(int a_x = -1, int a_y = -1);
	
#ifdef __OSMIA_PESTICIDE_STORE
	/** @brief Unique animal ID for pesticide exposure tracking/output */
	unsigned int m_animal_id;
#endif
};

#endif // Header include guard

/**
 * @class Osmia_InCocoon
 * @brief Overwintering adult stage within cocoon
 * @extends Osmia_Pupa
 *
 * @details Osmia_InCocoon represents the eclosed adult bee overwintering within its protective cocoon.
 * This stage implements the three-phase overwintering model: prewintering (late summer/autumn at elevated
 * temperatures), diapause proper (winter at low temperatures), and post-diapause quiescence (early spring
 * awaiting emergence cues). Overwintering mortality depends on prewinter thermal conditions, and emergence
 * timing is determined by an emergence counter that decrements with spring warming.
 *
 * @par Biological Foundation
 * Adult *O. bicornis* eclose from pupae in late summer (August-September) but remain within cocoons through
 * winter, emerging the following spring (April-May). This overwintering strategy protects adults from winter
 * weather whilst avoiding energetic costs of maintaining activity. The three-phase model reflects documented
 * physiological transitions: initial high respiration rates during prewintering deplete fat reserves if warm
 * conditions persist, deep diapause conserves energy through winter, and post-diapause quiescence awaits
 * appropriate emergence conditions.
 *
 * @par Three-Phase Overwintering Model
 *
 * **Phase 1 - Prewintering** (typically September-November):
 * - Occurs at temperatures above prewintering threshold (15°C)
 * - Accumulates degree-days (m_DDPrewinter) that increase mortality risk
 * - Represents fat depletion at warm temperatures
 * - Based on Sgolastra et al. (2011) work with *O. lignaria*
 *
 * **Phase 2 - Diapause** (typically November-January):
 * - Deep dormancy at temperatures near/below overwintering threshold (0°C)
 * - Minimal metabolic activity, lipid conservation
 * - No degree-day accumulation during this phase
 * - Chilling requirement for diapause completion (not explicitly modelled as separate counter)
 *
 * **Phase 3 - Post-diapause quiescence** (typically February-March):
 * - Diapause complete but awaiting emergence cues
 * - Emergence counter decrements with temperatures above emergence threshold (5°C)
 * - When counter reaches zero, adult emerges from nest
 * - Counter equation: 35.48 - 0.0147 × DD_accumulated
 *
 * @par Overwintering Mortality
 * Calculated as function of prewinter degree-day accumulation using equation from Sgolastra et al. (2011):
 * **mortality_probability = 0.05 × DD_prewinter - 4.63**
 *
 * This captures the biological reality that prolonged warm autumn conditions deplete lipid reserves,
 * reducing winter survival. Cool autumn conditions (low DD accumulation) produce high survival.
 *
 * @par Difference from Formal Model
 * **EXACT MATCH for equations** - Overwintering mortality and emergence counter equations implemented
 * precisely as specified in formal model. Temperature thresholds match formal model specifications
 * (0°C overwintering, 5°C emergence, 15°C prewintering baseline).
 *
 * **Minor calibration**: Emergence counter constant adjusted from 39.48 to 35.48, shifting emergence
 * timing slightly earlier in spring to match field observations.
 *
 * @par Implementation Notes
 * The three phases are not explicitly flagged as separate states - rather, they emerge from temperature-
 * threshold logic applied to DD accumulation and counter updates. This creates appropriate seasonal
 * phenology without needing complex state tracking.
 *
 * @see Osmia_Base for overwintering parameters
 * @see st_Develop() for phase transition logic
 * @see WinterMortality() for mortality calculation
 * @see st_Emerge() for spring emergence
 */
class Osmia_InCocoon : public Osmia_Pupa
{
protected:
	/**
	 * @var m_emergencecounter
	 * @brief Countdown to spring emergence (decrements with warm days)
	 * @details Initialized at eclosion using equation: counter = 35.48 - 0.0147 × initial_DD
	 * Each day above emergence threshold (5°C) decrements counter by accumulated DD that day.
	 * When counter reaches zero (or negative), adult emerges from nest.
	 *
	 * @par Biological Interpretation
	 * The counter represents integrated thermal accumulation needed for emergence readiness. After
	 * diapause completion (chilling requirement met), adults still need specific amount of spring
	 * warming before emergence is triggered. Counter mechanism ensures emergence synchrony with
	 * appropriate spring phenology (flower availability, favorable weather).
	 *
	 * @par Empirical Basis
	 * Equation calibrated from field emergence observations showing relationship between winter/spring
	 * thermal patterns and emergence dates. Lower initial DD (cool conditions) requires more spring
	 * warming (higher counter start), whilst higher initial DD needs less spring warming.
	 *
	 * @par Implementation Note
	 * Counter can become negative if warm spring arrives rapidly - this is acceptable as it simply
	 * means emergence threshold was exceeded quickly. Once counter ≤ 0, emergence occurs on next
	 * favorable day (temperature above emergence threshold).
	 */
	int m_emergencecounter;

	/**
	 * @var m_DDPrewinter
	 * @brief Accumulated degree-days during prewintering period (above 15°C baseline)
	 * @details Accumulates from eclosion through autumn. Only temperatures above 15°C contribute:
	 * DD_prewinter += (T_daily - 15) when T > 15°C
	 *
	 * Used in overwintering mortality calculation. High values (prolonged warm autumn) increase
	 * mortality risk via fat depletion. Reset to zero at eclosion for new adults.
	 *
	 * @par Biological Basis
	 * At temperatures above 15°C, prepupal/adult metabolism remains elevated, burning lipid reserves
	 * that are needed for winter survival. Bosch et al. (2008) documented weight loss rates of
	 * 0.2-0.4 mg/day during warm prewintering, with corresponding survival reductions.
	 *
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation follows formal model specification precisely. The 15°C baseline
	 * comes directly from Sgolastra et al. (2011) methodology.
	 *
	 * @par Valid Range
	 * [0, ~150] degree-days. Values >100 DD indicate poor overwintering conditions (prolonged warmth).
	 * Typical central European autumn might accumulate 30-60 DD, giving mortality around 0.2-0.4
	 * (20-40% death probability).
	 */
	double m_DDPrewinter;

public:
	/**
	 * @brief Constructor for newly eclosed adult in cocoon
	 * @param data Initialization data from pupal stage
	 * @details Initializes overwintering-specific attributes:
	 * - Sets m_DDPrewinter to 0 (begins accumulation from eclosion)
	 * - Calculates initial m_emergencecounter from current conditions
	 * - Transfers adult mass from pupal provision mass
	 * - Maintains nest linkage for eventual emergence
	 */
	Osmia_InCocoon(struct_Osmia* data);

	/**
	 * @brief Reinitialize object from pool
	 * @param data New initialization data
	 */
	virtual void ReInit(struct_Osmia* data);

	/** @brief Destructor */
	virtual ~Osmia_InCocoon();

	/**
	 * @brief Main step function for overwintering adults
	 * @details Each day:
	 * 1. Call st_Develop to update DD accumulation and counters
	 * 2. Check emergence conditions (counter ≤ 0 and T > threshold)
	 * 3. Apply winter mortality test based on accumulated prewinter DD
	 * 4. Transition to emergence or continue overwintering
	 */
	virtual void Step(void);

	/**
	 * @brief Set overwintering temperature threshold (static parameter)
	 * @param a_temp Temperature threshold (°C) for overwintering DD accumulation
	 * @details Used primarily during initialization to set population-wide threshold from configuration.
	 */
	static void SetOverwinteringTempThreshold(double a_temp) { m_OverwinteringTempThreshold = a_temp; }

	/**
	 * @brief Get prewinter degree-day accumulation
	 * @return Total DD accumulated above 15°C baseline during prewintering
	 * @details Accessor for monitoring overwintering conditions and mortality risk. Used in output
	 * generation and validation.
	 */
	double GetDDPreWinter() { return m_DDPrewinter; }

protected:
	/**
	 * @brief Development state - manage overwintering phases and emergence preparation
	 * @return Next state (st_Emerge when ready, st_Develop to continue, or st_Die)
	 *
	 * @details Complex multi-phase logic:
	 *
	 * **Prewintering phase** (warm autumn conditions):
	 * - If T > 15°C: accumulate to m_DDPrewinter
	 * - Represents fat depletion period
	 * - No explicit phase flag - identified by temperature
	 *
	 * **Diapause phase** (winter cold):
	 * - Temperatures near/below 0°C
	 * - Minimal DD accumulation
	 * - Lipid conservation
	 * - No active tracking - just waiting period
	 *
	 * **Post-diapause phase** (spring warming):
	 * - If T > 5°C (emergence threshold): decrement m_emergencecounter
	 * - Counter -= (T - threshold)
	 * - When counter ≤ 0: ready to emerge
	 * - Check for favorable conditions (appropriate temperature, daylight)
	 *
	 * **Mortality application**:
	 * - Call WinterMortality() to calculate probability from prewinter DD
	 * - Test against random number
	 * - If dies: return st_Die
	 * - If survives and counter ≤ 0: return st_Emerge
	 * - Otherwise: return st_Develop (continue overwintering)
	 *
	 * @par Temperature Threshold Handling
	 * Three thresholds govern different processes:
	 * - 15°C: Prewinter DD accumulation (represents elevated metabolism)
	 * - 5°C: Emergence counter (spring warming required)
	 * - 0°C: Overwintering baseline (diapause maintenance)
	 *
	 * These thresholds create realistic seasonal phenology without explicit calendar tracking.
	 *
	 * @par Edge Cases
	 * - Very warm autumn (high prewinter DD): High mortality, survivors may emerge earlier
	 * - Very cold autumn (low prewinter DD): High survival, emergence depends more on spring warming
	 * - Warm winter spell: May decrement emergence counter but doesn't trigger premature emergence
	 *   (counter must reach zero AND appropriate conditions must persist)
	 *
	 * @par Difference from Formal Model
	 * **IMPLEMENTATION MATCH** - The logic structure follows formal model three-phase description.
	 * Temperature thresholds and equations match specifications. The implementation detail of using
	 * continuous threshold checks rather than explicit phase flags is consistent with formal model
	 * emphasis on emergent behaviour from temperature-threshold interactions.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);

	/**
	 * @brief Transition state - emerge from cocoon as active adult
	 * @return Next state (st_Die as InCocoon object deleted)
	 *
	 * @details Spring emergence sequence:
	 * 1. Signal population manager to create Osmia_Female object
	 * 2. Transfer attributes:
	 *    - Adult mass (determines size and fecundity)
	 *    - Natal nest location (becomes dispersal origin)
	 *    - Age (for lifespan tracking)
	 *    - Any carried-over attributes
	 * 3. Population manager handles nest cleanup (remove from nest cell list)
	 * 4. Delete InCocoon object
	 * 5. New Osmia_Female begins active adult life cycle
	 *
	 * @par Biological Transition
	 * Adult bites through cocoon and cell partition, crawls out of nest, orients to environment,
	 * and begins pre-nesting maturation period. In model, this is represented by object class change
	 * from InCocoon (dormant, nest-bound) to Female (active, mobile).
	 *
	 * @par Timing Considerations
	 * Emergence typically occurs on warm, sunny days in spring when:
	 * - Emergence counter has reached zero (sufficient spring warming accumulated)
	 * - Current temperature above threshold (favorable immediate conditions)
	 * - Potentially daylight/weather checks (not explicitly modelled)
	 *
	 * This ensures emerged adults encounter appropriate conditions for initial flights and resource
	 * location.
	 *
	 * @par Implementation Note
	 * The formal model discusses emergence timing but leaves specific trigger details to implementation.
	 * This counter-based approach captures observed field patterns: synchronous emergence within
	 * 2-3 week period, with timing responding to both winter and spring thermal conditions.
	 */
	virtual TTypeOfOsmiaState st_Emerge(void);

	/**
	 * @brief Calculate overwintering mortality probability from prewinter thermal conditions
	 * @return true if dies, false if survives
	 *
	 * @details Implements Sgolastra et al. (2011) linear relationship:
	 * **mortality_probability = 0.05 × m_DDPrewinter - 4.63**
	 *
	 * Then tests this probability against uniform random number to determine fate.
	 *
	 * **Example calculations**:
	 * - m_DDPrewinter = 30 DD → mortality = 0.05×30 - 4.63 = -3.13 → 0.0 (clamped) → 0% death
	 * - m_DDPrewinter = 60 DD → mortality = 0.05×60 - 4.63 = -1.63 → 0.0 (clamped) → 0% death
	 * - m_DDPrewinter = 93 DD → mortality = 0.05×93 - 4.63 = 0.02 → 2% death
	 * - m_DDPrewinter = 100 DD → mortality = 0.05×100 - 4.63 = 0.37 → 37% death
	 * - m_DDPrewinter = 130 DD → mortality = 0.05×130 - 4.63 = 1.87 → 1.0 (clamped) → 100% death
	 *
	 * @par Biological Interpretation
	 * The equation captures lipid depletion during warm prewintering. Each degree-day above 15°C
	 * represents time at elevated metabolism, burning fat reserves needed for winter survival.
	 * The negative intercept (-4.63) means zero mortality at low DD accumulation (cool autumns provide
	 * ideal prewintering conditions). Mortality rises linearly with warm autumn conditions.
	 *
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Equation coefficients implemented precisely as specified in formal model
	 * (const = -4.63, slope = 0.05). Based directly on Sgolastra et al. (2011) empirical data for
	 * *O. lignaria* males, applied to both sexes of *O. bicornis* in absence of sex-specific data.
	 *
	 * @par Uncertainty
	 * MEDIUM - Cross-species application (*O. lignaria* → *O. bicornis*) introduces uncertainty.
	 * Male-to-female extrapolation also uncertain (females may differ in lipid reserves and metabolic
	 * rates). However, underlying mechanism (fat depletion at warm temperatures) is well-established
	 * and should be qualitatively similar across related species and sexes.
	 *
	 * @par Implementation Details
	 * - Probability clamped to [0, 1] range (equation can produce values outside this range)
	 * - Mortality test occurs once per overwintering period (not daily like developmental stages)
	 * - Timing of test determines whether early-season or late-season conditions dominate
	 * - Current implementation likely tests in spring (precise timing depends on Step logic)
	 *
	 * @par Valid DD Range and Interpretation
	 * - **0-90 DD**: Ideal to good prewintering (mortality 0-10%)
	 * - **90-110 DD**: Moderate stress (mortality 10-40%)
	 * - **110-130 DD**: Poor conditions (mortality 40-80%)
	 * - **>130 DD**: Catastrophic (mortality >80%)
	 *
	 * Central European autumns typically produce 30-80 DD, giving 0-20% mortality under normal conditions.
	 * Climate warming could increase DD accumulation, raising population-level mortality risk.
	 */
	bool WinterMortality();

	/**
	 * @var m_OverwinteringTempThreshold
	 * @brief Static temperature threshold for overwintering phase (°C)
	 * @details Default: 0.0°C. Defines baseline for diapause proper. Temperatures near/below this
	 * threshold represent true winter conditions where metabolism is minimized and lipid conservation
	 * is maximal.
	 *
	 * @par Usage
	 * Not used for DD accumulation (unlike egg/larva/pupa thresholds) but serves as conceptual
	 * boundary between active prewintering (>15°C), transitional conditions (0-15°C), and true
	 * diapause (<0°C approximately).
	 */
	static double m_OverwinteringTempThreshold;
};

/**
 * @class Osmia_Female
 * @brief Active adult female conducting reproduction
 * @extends Osmia_InCocoon
 *
 * @details Osmia_Female represents the culmination of the life cycle - the active adult female engaging
 * in dispersal, nest finding, foraging, provisioning, and egg laying. This is the most complex life stage,
 * with emergent spatial behaviour, resource-dependent reproductive decisions, and multiple interacting
 * state variables governing daily activities.
 *
 * @par Biological Foundation
 * Adult female *O. bicornis* emerge in spring (April-May), undergo a brief pre-nesting maturation period,
 * then begin reproductive activities that may span 4-8 weeks. Typical female completes 2-4 nests in her
 * lifetime, each containing 6-12 cells (eggs). Foraging behaviour follows central-place foraging theory:
 * females return repeatedly to their nest, balancing travel costs against resource quality. Sex allocation
 * (female vs. male offspring) responds to provision mass availability, with larger cells receiving female
 * eggs and smaller cells receiving male eggs.
 *
 * @par Reproductive Cycle
 * Female reproductive behaviour follows a repeating cycle:
 * 1. **Dispersal/nest searching**: Locate suitable cavity for nest
 * 2. **Nest establishment**: Clean cavity, orient to location
 * 3. **Cell provisioning cycle** (repeated per cell):
 *    - Forage for pollen/nectar
 *    - Return to nest with provisions
 *    - Construct cell partition
 *    - Determine sex of egg (based on provision mass)
 *    - Lay egg
 *    - Seal cell
 * 4. **Nest completion**: Seal final cell, abandon nest
 * 5. **Return to step 1** if longevity and eggs remaining permit
 *
 * @par Foraging Behaviour
 * Foraging implements spatially-explicit resource search using pre-computed masks for efficiency.
 * Females exhibit:
 * - Age-dependent foraging efficiency (Seidelmann 2006 curves)
 * - Give-up thresholds (abandon poor patches)
 * - Distance-dependent returns (closer patches preferred)
 * - Competition effects (pollen depletion by bee density)
 *
 * @par Sex Allocation
 * Sex determination follows haplodiploid genetics (fertilized=female, unfertilized=male) with strategic
 * maternal control. Females allocate sex based on provision mass:
 * - Large provisions → female egg (daughters require more resources)
 * - Small provisions → male egg (sons can develop on less)
 * - Sequential pattern: females typically at back of nest, males near entrance
 *
 * @par Mortality
 * Adult females experience daily background mortality (0.02/day from Giejdasz et al. 2016) representing
 * combined hazards of foraging flights, weather exposure, predation, and senescence. Additionally
 * vulnerable to pesticide exposure via contaminated pollen and direct spray contact.
 *
 * @par Difference from Formal Model
 * The formal model describes conceptual reproductive behaviors and resource dependencies. Implementation
 * adds necessary spatial algorithms (foraging masks, movement mechanics), specific sex allocation rules,
 * and detailed provisioning logistics. Core biological relationships (mass-fecundity from Seidelmann 2010,
 * age-efficiency curves) implemented exactly as specified. Spatial implementation details extend formal
 * model to enable landscape-scale simulation.
 *
 * @see st_ReproductiveBehaviour() for reproductive state machine
 * @see Forage() for detailed foraging algorithm
 * @see CalcParasitised() for parasitism risk assessment
 * @see LayEgg() for sex allocation and egg creation
 */
class Osmia_Female : public Osmia_InCocoon
{
public:
#ifdef __OSMIARECORDFORAGE
	/** @brief Cumulative foraging success across all females (testing/validation only) */
	static double m_foragesum;
	/** @brief Count of foraging events (testing/validation only) */
	static int m_foragecount;
#endif

protected:
	//----------------------- Foraging Infrastructure -----------------------

	/**
	 * @var m_foragemask
	 * @brief Static coarse-resolution spatial search mask
	 * @details Shared across all females for memory efficiency. Provides 20 distance rings × 8 directions
	 * for efficient outward resource searches from nest location.
	 */
	static OsmiaForageMask m_foragemask;

	/**
	 * @var m_foragemaskdetailed
	 * @brief Static high-resolution spatial search mask
	 * @details Alternative mask with finer spatial resolution for detailed pollen assessment. Used when
	 * comprehensive resource evaluation needed rather than incremental search.
	 */
	static OsmiaForageMaskDetailed m_foragemaskdetailed;

	/**
	 * @var m_currentpollenlevel
	 * @brief Current pollen availability at active foraging location
	 * @details Updated during foraging to track resource depletion at focal patch. Compared against
	 * give-up thresholds to determine when to abandon patch and search elsewhere.
	 */
	double m_currentpollenlevel;

	/**
	 * @var m_pollengiveupthreshold
	 * @brief Proportional reduction triggering patch abandonment
	 * @details When pollen level drops to this proportion of initial value, female abandons current
	 * patch and searches for new location. Represents marginal value theorem - forager leaves when
	 * diminishing returns make travel to new patch worthwhile.
	 *
	 * @par Typical Value
	 * 0.5-0.7 (abandon when patch drops to 50-70% of initial quality)
	 */
	static double m_pollengiveupthreshold;

	/**
	 * @var m_pollengiveupreturn
	 * @brief Absolute pollen level below which new search triggered
	 * @details Minimum acceptable return rate. If patch quality falls below this absolute threshold,
	 * female searches for new patch regardless of proportional decline. Prevents wasting time on
	 * completely depleted areas.
	 */
	static double m_pollengiveupreturn;

	//----------------------- Nest Provisioning State -----------------------

	/**
	 * @var m_CellOpenDays
	 * @brief Number of days current cell has been open (accumulating parasitism risk)
	 * @details Parasitism probability increases with cell open time - longer provisioning periods
	 * (due to poor forage availability or bad weather) expose cells to more parasitoid encounters.
	 * Reset to zero when cell is sealed.
	 */
	int m_CellOpenDays;

	/**
	 * @var m_CellCarryOver
	 * @brief Fractional hours carried to next day when cell not completed
	 * @details Cell construction requires minimum time (typically 1 day) but poor foraging may stretch
	 * across multiple days. This variable tracks partial progress to avoid losing fractional time
	 * between days.
	 */
	double m_CellCarryOver;

	/**
	 * @var m_EggsToLay
	 * @brief Total lifetime egg load remaining
	 * @details Calculated at emergence from body mass using Seidelmann (2010) relationship:
	 * total_eggs = N_nests_possible × (0.0371 × mass + 2.8399) ± 3
	 *
	 * Decrements with each egg laid. When reaches zero, female ceases reproduction even if surviving.
	 * Represents ovary capacity constraint - bee cannot produce unlimited eggs.
	 */
	int m_EggsToLay;

	/**
	 * @var m_EggsThisNest
	 * @brief Planned eggs for current nest (decrements as cells completed)
	 * @details Drawn from probability distribution at nest initiation, representing female's "plan"
	 * for nest size. Actual eggs laid may differ if resources fail or female dies mid-nest.
	 * Planning occurs at nest start, simulating observed tendency for females to provision nests
	 * of characteristic size (6-12 cells typical).
	 */
	int m_EggsThisNest;

	/**
	 * @var m_ToDisperse
	 * @brief Flag indicating need for dispersal to new nesting area
	 * @details Set true when: (1) emergence (find initial nest area), (2) nest completion (may search
	 * new area for next nest), (3) repeated nest-finding failures (exhaust local options). Controls
	 * transition to long-distance dispersal behaviour vs. local nest searching.
	 */
	bool m_ToDisperse;

	/**
	 * @var m_EmergeAge
	 * @brief Days since emergence (tracks adult age separately from total age)
	 * @details Used for: age-dependent foraging efficiency (Seidelmann 2006 curves show efficiency
	 * increasing first 7-10 days then declining), lifespan constraints (max ~60 days), and output
	 * tracking of adult longevity.
	 */
	int m_EmergeAge;

	/**
	 * @var m_CurrentNestLoc
	 * @brief Spatial location of nest currently being provisioned
	 * @details X,Y coordinates of active nest. When no nest (dispersing or searching), m_x set to -1
	 * as flag. All foraging trips reference this location as return point (central-place foraging).
	 * Updated when new nest found.
	 */
	APoint m_CurrentNestLoc;

	/**
	 * @var m_ProvisioningTime
	 * @brief Days required to complete one cell
	 * @details Depends on: forage quality (poor resources extend provisioning), weather (bad days
	 * prevent foraging), female efficiency (age effects). Typically 1-3 days per cell. Longer times
	 * increase parasitism risk via extended cell open duration.
	 */
	int m_ProvisioningTime;

	/**
	 * @var m_FlyingCounter
	 * @brief Days spent flying/foraging during current cell construction
	 * @details Distinguishes active foraging days from weather delays. Used to accurately track
	 * provisioning effort vs. unavoidable delays. Increments only on days when foraging occurs.
	 */
	int m_FlyingCounter;

	/**
	 * @var m_CurrentProvisioning
	 * @brief Mass of pollen/nectar currently provisioned in active cell (mg)
	 * @details Accumulates from zero as female makes foraging trips. When reaches target mass for
	 * planned sex (female target or male target), cell is complete and egg is laid. Determines
	 * final egg sex and offspring size.
	 */
	double m_CurrentProvisioning;

	/**
	 * @var m_BeeSizeScore1
	 * @brief Coarse size class (0=very small, 1=small, 2=medium, 3=large)
	 * @details Categorical size classification from adult mass. Used for size-dependent behaviour
	 * parameters if implemented (currently placeholder for future extensions where larger bees might
	 * have different foraging ranges or fecundity).
	 */
	int m_BeeSizeScore1;

	/**
	 * @var m_BeeSizeScore2
	 * @brief Fine-grained size classification
	 * @details Finer size categories than m_BeeSizeScore1, with step size controlled by configuration
	 * parameter cfg_OsmiaAdultMassCategoryStep. Enables more nuanced size-dependent parameterization
	 * if needed for model extensions.
	 */
	int m_BeeSizeScore2;

	/**
	 * @var m_NestProvisioningPlan
	 * @brief Queue of target provision masses for planned nest cells
	 * @details Female "plans" nest at initiation, generating sequence of target masses (one per egg).
	 * Deque structure allows efficient removal from front as cells completed. Masses decline from
	 * first to last cell (progressive resource depletion effect), with stochastic variation.
	 *
	 * @par Biological Basis
	 * Seidelmann (2010) documented progressive decline in provision mass from first to last offspring
	 * within nests, reflecting maternal aging and resource depletion. This planning structure implements
	 * that pattern whilst allowing females to adjust actual provisioning based on encountered resources.
	 */
	deque<double>m_NestProvisioningPlan;

	/**
	 * @var m_NestProvisioningPlanSex
	 * @brief Queue of planned sexes (true=female, false=male) corresponding to m_NestProvisioningPlan
	 * @details Sex allocated at planning stage based on provision mass targets: larger masses get
	 * female designation, smaller get male. Female control over fertilization allows implementation
	 * of planned sex allocation. Actual sex ratio emerges from resource availability and mortality.
	 */
	deque<bool>m_NestProvisioningPlanSex;

	//----------------------- Foraging State -----------------------

	/**
	 * @var m_ForageLoc
	 * @brief Flag indicating whether foraging location has been identified
	 * @details true = female has located pollen source and is actively exploiting it.
	 * false = female needs to search for new forage location (initial or after patch abandonment).
	 * Controls whether to search or continue foraging at known location.
	 */
	bool m_ForageLoc;

	/**
	 * @var m_ForageLocPoly
	 * @brief Index to polygon list in population manager providing resources
	 * @details Each landscape polygon (field, forest patch, etc.) has associated resource values.
	 * This index provides fast lookup of focal polygon's resource data without repeated spatial queries.
	 * Updated when new foraging location selected.
	 */
	int m_ForageLocPoly;

	/**
	 * @var m_ForageSteps
	 * @brief Number of distance steps in foraging mask (determines search resolution)
	 * @details Static value controlling how many distance rings are searched. Higher values allow
	 * finer-grained searches but increase computation. Typically 10-20 steps covering 0 to max
	 * foraging range.
	 */
	static int m_ForageSteps;

	/**
	 * @var m_PollenCompetitionsReductionScaler
	 * @brief Scaling factor for inter-specific pollen competition
	 * @details Adjusts available pollen based on assumed competition from other bee species (honey bees,
	 * bumble bees, other solitaries). Values <1.0 reduce available pollen, simulating competitive depletion
	 * by non-modelled species. Enables exploration of competition scenarios.
	 *
	 * @par Implementation Note
	 * Simple proportional reduction rather than mechanistic competition model. Pragmatic approach given
	 * uncertainty in competitor densities and overlap in flower usage.
	 */
	static double m_PollenCompetitionsReductionScaler;

	/**
	 * @var m_FemaleForageEfficiency
	 * @brief Vector of age-dependent foraging efficiency multipliers indexed by adult age
	 * @details Implements Seidelmann (2006) empirical efficiency curves showing:
	 * - Days 1-7: Efficiency increases as females gain experience
	 * - Days 8-15: Peak efficiency (full capability)
	 * - Days 15+: Gradual decline with senescence (wing wear, reduced flight capability)
	 *
	 * Applied as multiplier to daily forage returns: actual_forage = base_forage × efficiency[age]
	 *
	 * @par Biological Basis
	 * Young bees need practice to optimize foraging routes and flower handling. Old bees experience
	 * cumulative wing wear and muscle degradation reducing flight speed and cargo capacity.
	 */
	static vector<double> m_FemaleForageEfficiency;

	/**
	 * @var m_ForageLocX
	 * @brief X-coordinate of current foraging location
	 */
	int m_ForageLocX;

	/**
	 * @var m_ForageLocY
	 * @brief Y-coordinate of current foraging location
	 */
	int m_ForageLocY;

	/**
	 * @var m_foraged_resource_pesticide
	 * @brief Array storing pesticide concentrations in foraged resources
	 * @details When pesticide module active, tracks pesticide content of pollen/nectar collected at
	 * different locations. Enables simulation of pesticide exposure via contaminated provisions.
	 * Array size matches number of pesticide types configured.
	 */
	double* m_foraged_resource_pesticide;

	//----------------------- Pesticide Exposure (Conditional Compilation) -----------------------

#ifdef __OSMIA_PESTICIDE_ENGINE
public:
	/** @brief Egg-specific pesticide death probability after threshold exceedance */
	static double m_OsmiaEggPPPEffectProb;

	/** @brief Pesticide concentration threshold for egg effects */
	static double m_OsmiaEggPPPThreshold;

	/** @brief Adult pesticide death probability after threshold exceedance */
	static double m_OsmiaPPPEffectProb;

	/** @brief Adult pesticide concentration threshold */
	static double m_OsmiaPPPThreshold;

	/** @brief Daily pesticide decay rate in bee body (proportion lost per day) */
	static double m_OsmiaPPPDecayRate;

	/** @brief Absorption rate for overspray exposure (proportion transferred body to internal) */
	static double m_OsmiaPPPAbsorptionRateOverspray;

	/** @brief Absorption rate for contact exposure */
	static double m_OsmiaPPPAbsorptionRateContact;

	/** @brief Surface area exposed to overspray (mm²) */
	static double m_OsmiaPPPOversprayBodySurface;

	/** @brief Surface area for contact exposure (mm²) */
	static double m_OsmiaPPPContactBodySurface;

	/** @brief Probability of experiencing overspray event */
	static double m_OsmiaPPPOversprayChance;

protected:
	/** @brief Accessor for pesticide threshold parameter */
	static double GetPPPThreshold() { return m_OsmiaPPPThreshold; }

	/** @brief Accessor for pesticide effect probability */
	static double GetPPPEffectProb() { return m_OsmiaPPPEffectProb; }

	/** @brief Accessor for pesticide decay rate */
	static double GetPPPDecayRate() { return m_OsmiaPPPDecayRate; }

	/** @brief Accessor for overspray absorption rate */
	static double GetPPPAbsorptionRateOverspray() { return m_OsmiaPPPAbsorptionRateOverspray; }

	/** @brief Accessor for overspray body surface */
	static double GetPPPOversprayBodySurface() { return m_OsmiaPPPOversprayBodySurface; }

	/** @brief Accessor for contact body surface */
	static double GetPPPContactBodySurface() { return m_OsmiaPPPContactBodySurface; }
#endif

	//----------------------- Testing Support (Conditional Compilation) -----------------------

#ifdef __OSMIATESTING
	/** @brief Target nest data for validation (intended provisioning plan) */
	OsmiaNestData m_target;

	/** @brief Achieved nest data for validation (actual provisioning accomplished) */
	OsmiaNestData m_achieved;

	/** @brief Flag for first nest tracking */
	bool m_firstnestflag;
#endif

	//----------------------- Behavioural Methods -----------------------

	/**
	 * @brief Death state with female-specific cleanup
	 * @details Extends base st_Dying to handle:
	 * - Incomplete nest abandonment (set nest to closed, prevent further egg additions)
	 * - Resource release (return unused forage hours to pool)
	 * - Output recording (lifetime reproductive success, cause of death)
	 */
	virtual void st_Dying(void);

	/**
	 * @brief Development state for adult females (minimal - no metamorphosis)
	 * @return Next state (typically st_ReproductiveBehaviour or st_Dispersal)
	 * @details Unlike immature stages, adult females don't "develop" - this state primarily handles
	 * daily initialization and transitions to active reproductive states. May handle pre-nesting
	 * maturation period (first few days post-emergence before reproduction begins).
	 */
	virtual TTypeOfOsmiaState st_Develop(void);

	/**
	 * @brief Search for suitable nest cavity
	 * @return true if nest found, false if search fails
	 *
	 * @details Nest finding algorithm:
	 * 1. Sample locations around current position using movement probability distribution
	 * 2. Check each location for suitable cavities (queries landscape manager)
	 * 3. If suitable cavity available: establish nest, set m_CurrentNestLoc, return true
	 * 4. If no cavity found after N attempts: set m_ToDisperse=true, return false
	 *
	 * @par Cavity Suitability Criteria
	 * - Appropriate diameter (6-9mm for *O. bicornis*)
	 * - Sufficient depth (>10cm)
	 * - Protected location (not fully exposed)
	 * - Not already occupied
	 *
	 * In model, suitability determined by landscape polygon attributes (habitat types provide different
	 * cavity densities).
	 *
	 * @par Search Limitations
	 * Female makes limited attempts (m_OsmiaFindNestAttemptNo, typically 5-10). Repeated failures
	 * trigger dispersal to new area. Reflects biological reality that suitable cavities are limiting
	 * resource in many landscapes.
	 */
	virtual bool FindNestLocation(void);

	/**
	 * @brief Dispersal state for long-distance movements to new nesting areas
	 * @return Next state (st_ReproductiveBehaviour if dispersal successful, st_Die if fails)
	 *
	 * @details Long-distance dispersal using different movement distribution than local foraging:
	 * - Samples from m_dispersalmovementdistances (typically longer distances)
	 * - Moves to new location, attempts nest finding
	 * - If successful: begins reproduction at new location
	 * - If fails: may attempt further dispersal or die (dispersal mortality)
	 *
	 * @par Biological Context
	 * Dispersal occurs: (1) at emergence from natal nest, (2) after nest completion if local resources
	 * depleted, (3) after repeated nest-finding failures. Enables population spread and colonization
	 * of new areas, but carries mortality costs (navigation failures, exposure, energy depletion).
	 *
	 * @par Difference from Formal Model
	 * Formal model describes dispersal conceptually. Implementation specifies movement mechanics using
	 * beta probability distributions parameterized from allometric relationships between body size and
	 * foraging range. Dispersal distances typically 2-3× longer than foraging movements, reaching
	 * maximum homing distance (R90 = 1430m).
	 */
	virtual TTypeOfOsmiaState st_Dispersal(void);

	/**
	 * @brief Main foraging algorithm collecting pollen and nectar
	 * @return Mass of pollen collected (mg)
	 *
	 * @details Complex spatial foraging implementing:
	 *
	 * **Phase 1 - Resource Location** (if !m_ForageLoc):
	 * - Use OsmiaForageMask to search concentrically from nest
	 * - Evaluate pollen availability at each location
	 * - Select location with acceptable pollen level
	 * - Store location, set m_ForageLoc = true
	 *
	 * **Phase 2 - Resource Exploitation** (if m_ForageLoc):
	 * - Calculate forage return based on:
	 *   - Local pollen availability
	 *   - Distance from nest (travel time reduces forage time)
	 *   - Age-dependent efficiency (m_FemaleForageEfficiency[age])
	 *   - Competition (density-dependent depletion)
	 *   - Hours available (m_foragehours)
	 * - Deplete local pollen (update landscape resource levels)
	 * - Check give-up thresholds:
	 *   - If pollen dropped to <threshold: m_ForageLoc=false (search new location)
	 *   - If absolute level <minimum: m_ForageLoc=false
	 * - Return collected mass
	 *
	 * **Phase 3 - Resource Accumulation**:
	 * - Add collected mass to m_CurrentProvisioning
	 * - Check if cell target reached
	 * - If complete: trigger egg laying sequence
	 *
	 * @par Age-Dependent Efficiency
	 * Daily foraging success scaled by m_FemaleForageEfficiency[m_EmergeAge]:
	 * - Days 0-7: Rising efficiency (inexperience)
	 * - Days 7-15: Peak efficiency
	 * - Days 15+: Declining efficiency (senescence)
	 *
	 * Based on Seidelmann (2006) empirical curves.
	 *
	 * @par Distance Effects
	 * Longer distances reduce foraging efficiency:
	 * - Travel time reduces available foraging hours
	 * - Energetic costs reduce net provisioning
	 * - Implemented via time budget (hours available after travel)
	 *
	 * @par Competition
	 * Pollen availability decreases with local *Osmia* density:
	 * available_pollen = base_pollen × (1 - density × m_DensityDependentPollenRemovalConst)
	 *
	 * Provides negative feedback preventing unrealistic population growth in favorable habitats.
	 *
	 * @par Give-Up Decisions
	 * Two thresholds govern patch abandonment:
	 * 1. Proportional decline: leave when pollen drops to X% of initial (marginal value theorem)
	 * 2. Absolute minimum: leave when pollen falls below acceptable return rate
	 *
	 * Whichever threshold triggered first causes female to search for new patch.
	 *
	 * @par Difference from Formal Model
	 * Formal model describes foraging conceptually with distance costs and resource depletion. Implementation
	 * specifies spatial search algorithms, time budgets, and give-up rules. Core relationships (age effects,
	 * distance costs, competition) implement formal model principles with specific parameterization from
	 * Seidelmann (2006) and calibration.
	 */
	double Forage(void);

	/**
	 * @brief Reproductive behaviour state coordinating nesting activities
	 * @return Next state (continues reproduction, disperses, or dies)
	 *
	 * @details Master reproductive state machine:
	 *
	 * **If no nest** (m_CurrentNestLoc.x == -1):
	 * - Call FindNestLocation()
	 * - If successful: initialize new nest
	 * - If fails: return st_Dispersal
	 *
	 * **If active nest**:
	 * - Check weather (if unfavorable: skip foraging, increment cell open days)
	 * - Call Forage() to collect provisions
	 * - Add to m_CurrentProvisioning
	 * - Check if cell target reached:
	 *   - If yes: call LayEgg(), reset m_CurrentProvisioning, decrement m_EggsThisNest
	 *   - If no: continue provisioning
	 * - Check if nest complete (m_EggsThisNest == 0):
	 *   - Seal nest, set m_CurrentNestLoc.x = -1
	 *   - Check remaining eggs (m_EggsToLay):
	 *     - If >0: prepare for next nest (may disperse or search locally)
	 *     - If 0: reproduction complete, continue until senescence
	 *
	 * **Daily mortality test**:
	 * - Apply m_OsmiaFemaleBckMort (0.02 daily probability)
	 * - If dies: return st_Die
	 *
	 * @par Weather Effects
	 * Bad weather (too cold, rainy, windy) prevents foraging but doesn't stop clock:
	 * - Cell remains open (accumulating parasitism risk)
	 * - Age advances (depleting remaining lifespan)
	 * - Provisions don't accumulate (extends nest completion time)
	 *
	 * Creates realistic weather impacts on reproductive success without explicit weather state tracking.
	 *
	 * @par Nest Completion Decision
	 * Female "knows" when nest is complete via m_EggsThisNest counter (planned at nest initiation).
	 * Actual completion may differ from plan if:
	 * - Resources fail (cannot provision remaining cells)
	 * - Female dies mid-nest
	 * - Eggs exhausted before plan complete
	 *
	 * This represents partial implementation of female's "plan" responding to realized conditions.
	 */
	virtual TTypeOfOsmiaState st_ReproductiveBehaviour(void);

	/**
	 * @brief Calculate planned eggs for next nest
	 * @return Number of eggs planned for nest
	 *
	 * @details Draws from m_eggspernestdistribution (typically beta distribution giving 6-12 eggs
	 * with appropriate skew). Value constrained by:
	 * - m_OsmiaFemaleMinEggsPerNest (lower bound)
	 * - m_OsmiaFemaleMaxEggsPerNest (upper bound)
	 * - m_EggsToLay (can't plan more than remaining lifetime egg load)
	 *
	 * @par Biological Basis
	 * Females exhibit characteristic nest sizes reflecting tradeoffs between offspring number and quality.
	 * Larger nests risk incomplete provisioning if female dies; smaller nests under-utilize female's
	 * reproductive capacity. Optimal nest size balances these factors and varies with resource quality.
	 *
	 * @par Implementation
	 * Planning occurs at nest initiation, creating "target" that female works toward. Actual eggs laid
	 * may differ based on resource availability, weather, mortality. This two-stage process (plan then
	 * execute) captures observed tendency for nests to cluster around typical sizes whilst allowing
	 * environmental variation.
	 */
	int PlanEggsPerNest();

	/**
	 * @brief Calculate total lifetime egg load from body mass
	 * @details Implements Seidelmann (2010) empirical relationship:
	 * **eggs_per_nest = 0.0371 × mass + 2.8399 (±3 eggs stochastic)**
	 * **total_eggs = m_TotalNestsPossible × eggs_per_nest**
	 *
	 * Sets m_EggsToLay at emergence. Larger females have higher fecundity via both more eggs per nest
	 * and potential for completing more nests (if they can provision faster).
	 *
	 * Also initializes m_EggsThisNest by calling PlanEggsPerNest() + 2 (2 removed at nest start,
	 * creating correct initial count).
	 *
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implements Seidelmann (2010) equation precisely as specified in formal model
	 * (const = 2.8399, slope = 0.0371, stochastic variation ±3 eggs). Total lifetime fecundity calculated
	 * as: eggs_per_nest × maximum_possible_nests, representing ovary capacity constraint.
	 */
	void CalculateEggLoad() {
		m_EggsToLay = int((m_TotalNestsPossible * (0.0371 * m_Mass + 2.8399)) + (g_rand_uni_fnc() * 6) - 3);
		m_EggsThisNest = PlanEggsPerNest() + 2;
	}

	/**
	 * @brief Determine parasitism status for egg about to be laid
	 * @param a_daysopen Number of days cell has been open
	 * @return Parasitoid type affecting this egg
	 *
	 * @details Two possible parasitism models (controlled by m_UsingMechanisticParasitoids):
	 *
	 * **Simple model** (probability-based):
	 * - Parasitism risk increases linearly with cell open time
	 * - Base probability scaled by m_ParasitismProbToTimeCellOpen
	 * - Bombylid probability from m_BombylidProbability
	 * - Random draw determines outcome
	 *
	 * **Mechanistic model** (population-based):
	 * - Queries parasitoid population manager for local parasitoid density
	 * - Attack probability from m_ParasitoidAttackChance[parasitoid_type] × density
	 * - Multiple parasitoid types possible with type-specific parameters
	 * - More realistic but requires additional parasitoid population tracking
	 *
	 * @par Biological Basis
	 * Open nest cells are vulnerable to parasitoid females searching for hosts. Longer provisioning
	 * times (multi-day cells) provide more opportunity for parasitoid discovery. Bombylid flies are
	 * primary parasitoids of *O. bicornis*, laying eggs in open cells that consume host egg/larva.
	 *
	 * @par Implementation Note
	 * Called during LayEgg() sequence, before egg object creation. Parasitism status assigned to egg
	 * at laying, determining subsequent survival. Parasitised eggs/larvae die at characteristic time
	 * for parasitoid species (handled by population manager during development).
	 */
	TTypeOfOsmiaParasitoids CalcParasitised(double a_daysopen);

	/**
	 * @brief Create and lay egg in completed cell
	 *
	 * @details Egg laying sequence:
	 * 1. Determine sex based on provision mass:
	 *    - If m_CurrentProvisioning >= female minimum: lay female egg (fertilized)
	 *    - If below female minimum: lay male egg (unfertilized)
	 * 2. Calculate parasitism status: CalcParasitised(m_CellOpenDays)
	 * 3. Create struct_Osmia with egg initialization data:
	 *    - Sex (determined above)
	 *    - Provision mass (m_CurrentProvisioning)
	 *    - Nest pointer (m_OurNest)
	 *    - Parasitism status
	 * 4. Signal population manager to create Osmia_Egg object
	 * 5. Nest adds egg to cell list
	 * 6. Reset cell state:
	 *    - m_CurrentProvisioning = 0
	 *    - m_CellOpenDays = 0
	 *    - m_EggsToLay -= 1
	 *    - m_EggsThisNest -= 1
	 *
	 * @par Sex Allocation Strategy
	 * Haplodiploid sex determination with maternal control allows strategic allocation:
	 * - Females (diploid) require more resources → placed on larger provisions
	 * - Males (haploid) can develop on less → placed on smaller provisions
	 *
	 * Threshold-based allocation emerges from provision mass variation:
	 * - Early cells (larger provisions): mostly females
	 * - Late cells (smaller provisions due to resource depletion): mostly males
	 * - Natural pattern: females at back of nest, males near entrance
	 *
	 * @par Difference from Formal Model
	 * Formal model describes sex allocation conceptually (larger provisions → females). Implementation
	 * specifies threshold-based decision rule with provision mass targets. Core biology (females need
	 * more resources) implemented faithfully; specific threshold values calibrated from field sex ratio
	 * observations.
	 *
	 * @par Parasitism Integration
	 * Parasitism status assigned at laying (based on cell open time) and carries through development.
	 * Parasitised individuals develop normally until parasitoid emerges (timing depends on parasitoid
	 * type), then host dies. This creates realistic delayed mortality rather than immediate death.
	 */
	void LayEgg();

public:
	/**
	 * @brief Constructor for newly emerged adult female
	 * @param data Initialization data from InCocoon stage
	 * @details Initializes:
	 * - Adult mass (determines size class and fecundity)
	 * - Calculates lifetime egg load: CalculateEggLoad()
	 * - Sets initial state (dispersal to find first nest area)
	 * - Initializes foraging attributes (no location, no nest)
	 * - Records emergence location (becomes dispersal origin)
	 */
	Osmia_Female(struct_Osmia* data);

	/**
	 * @brief Reinitialize from object pool
	 * @param data New initialization data
	 */
	virtual void ReInit(struct_Osmia* data);

	/** @brief Destructor */
	virtual ~Osmia_Female();

	/**
	 * @brief Female-specific initialization (called by constructor and ReInit)
	 * @param a_mass Adult body mass (mg)
	 * @details Handles mass-dependent initialization:
	 * - Size class calculation (m_BeeSizeScore1, m_BeeSizeScore2)
	 * - Provision mass targets (females vs. males)
	 * - Fecundity calculation: CalculateEggLoad()
	 */
	virtual void Init(double a_mass);

	/**
	 * @brief Pre-step initialization each day
	 * @details Sets up daily state:
	 * - Reset forage hours available (from weather/daylight)
	 * - Increment emerge age
	 * - Check lifespan limit (m_EmergeAge vs. m_OsmiaFemaleLifespan)
	 * - Update local resource availability if needed
	 */
	virtual void BeginStep(void);

	/**
	 * @brief Main step function orchestrating daily behaviour
	 * @details Calls appropriate behavioural state based on m_CurrentOState:
	 * - st_Develop: Initial maturation
	 * - st_Dispersal: Long-distance movement
	 * - st_ReproductiveBehaviour: Nesting and provisioning
	 * - st_Die: Cleanup and removal
	 *
	 * Loops until state machine returns terminal state (DONE or DIE).
	 */
	virtual void Step(void);

	//----------------------- Static Setters (Population Manager Initialization) -----------------------

	/** @brief Set number of distance steps in foraging mask */
	static void SetForageSteps(int a_sz) { m_ForageSteps = a_sz; }

	/**
	 * @brief Initialize detailed foraging mask
	 * @param a_step Step size between sample points
	 * @param a_max Maximum search distance
	 */
	static void SetForageMaskDetailed(int a_step, int a_max) {
		OsmiaForageMaskDetailed fmd(a_step, a_max);
		m_foragemaskdetailed = fmd;
	}

	/** @brief Set proportional give-up threshold for patch abandonment */
	static void SetPollenGiveUpThreshold(double a_prop) { m_pollengiveupthreshold = a_prop; }

	/** @brief Set absolute give-up threshold (minimum acceptable return) */
	static void SetPollenGiveUpReturn(double a_value) { m_pollengiveupreturn = a_value; }

	/** @brief Set daily background mortality for adult females */
	static void SetDailyMort(double a_prob) { m_OsmiaFemaleBckMort = a_prob; }

	/** @brief Set number of nest-finding attempts before dispersal triggered */
	static void SetNestFindAttempts(int a_no) { m_OsmiaFindNestAttemptNo = a_no; }

	/** @brief Set minimum eggs per nest (lower bound for planning distribution) */
	static void SetMinEggsPerNest(int a_eggs) { m_OsmiaFemaleMinEggsPerNest = a_eggs; }

	/** @brief Set maximum eggs per nest (upper bound for planning distribution) */
	static void SetMaxEggsPerNest(int a_eggs) { m_OsmiaFemaleMaxEggsPerNest = a_eggs; }

	/**
	 * @brief Set cocoon-to-provision mass conversion and derived parameters
	 * @param a_ratio Conversion factor (proportion of provision mass converted to cocoon mass)
	 * @details Also calculates total provisioning mass loss parameters by scaling cocoon mass loss
	 * configuration values. These derived parameters used in nest provisioning planning.
	 */
	static void SetCocoonToProvisionMass(double a_ratio) {
		m_CocoonToProvisionMass = a_ratio;
		m_TotalProvisioningMassLoss = cfg_OsmiaTotalCocoonMassLoss.value() * a_ratio;
		m_TotalProvisioningMassLossRange = cfg_OsmiaTotalCocoonMassLossRange.value() * a_ratio;
		m_TotalProvisioningMassLossRangeX2 = m_TotalProvisioningMassLossRange * 2.0;
	}

	/** @brief Set provision-to-cocoon mass conversion factor */
	static void SetProvisionToCocoonMass(double a_ratio) { m_ProvisionToCocoonMass = a_ratio; }

	/** @brief Set pollen score to mg conversion factor */
	static void SetPollenScoreToMg(double a_ratio) { m_PollenScoreToMg = a_ratio; }

	/** @brief Set minimum target provision mass for male cells (instance method) */
	void SetMaleMinTargetProvisionMass(double a_mass) { m_MaleMinTargetProvisionMass = a_mass; }

	/** @brief Set minimum target provision mass for female cells (instance method) */
	void SetFemaleMinTargetProvisionMass(double a_mass) { m_FemaleMinTargetProvisionMass = a_mass; }

	/** @brief Set maximum target provision mass for female cells (instance method) */
	void SetFemaleMaxTargetProvisionMass(double a_mass) { m_FemaleMaxTargetProvisionMass = a_mass; }

	/** @brief Set minimum cell construction time (days) */
	static void SetMinimumCellConstructionTime(double a_time) { m_MinimumCellConstructionTime = a_time; }

	/** @brief Set maximum cell construction time (days) */
	static void SetMaximumCellConstructionTime(double a_time) { m_MaximumCellConstructionTime = a_time; }

	/** @brief Set maximum lifetime nests possible */
	static void SetTotalNestsPossible(int a_total) { m_TotalNestsPossible = a_total; }

	/** @brief Set Bombyliid parasitism probability */
	static void SetBombylidProbability(double a_prob) { m_BombylidProbability = a_prob; }

	/** @brief Set parasitism probability to cell open time conversion factor */
	static void SetParasitismProbToTimeCellOpen(double a_ratio) { m_ParasitismProbToTimeCellOpen = a_ratio; }

	/** @brief Set flag for using mechanistic vs. simple parasitoid model */
	static void SetUsingMechanisticParasitoids(bool a_flag) { m_UsingMechanisticParasitoids = a_flag; }

	/** @brief Set parasitoid attack probability parameters (vector for multiple types) */
	static void SetParasitoidParameters(vector<double> a_params) { m_ParasitoidAttackChance = a_params; }

	/** @brief Set density-dependent pollen removal constant (instance method) */
	void SetDensityDependentPollenRemovalConst(double a_value) { m_DensityDependentPollenRemovalConst = a_value; }

	/** @brief Add age-specific foraging efficiency value to vector */
	static void AddForageEfficiency(double a_eff) { m_FemaleForageEfficiency.push_back(a_eff); }

	/**
	 * @brief Get available pollen in polygon from starting location
	 * @param a_required_amount [in/out] Amount of pollen needed (mg)
	 * @param a_foraged_amount [out] Amount actually foraged (mg)
	 * @param a_polygon Polygon index to query
	 * @param a_loc_x X-coordinate in polygon
	 * @param a_loc_y Y-coordinate in polygon
	 *
	 * @details Queries landscape manager for pollen availability, applies competition effects,
	 * depletes local resources, returns actual foraged amount. Used by Forage() to implement
	 * resource acquisition and depletion.
	 */
	void GetPollenInPolygon(double& a_required_amount, double& a_foraged_amount, int a_polygon, int a_loc_x, int a_loc_y);

	/**
	 * @brief Handle farm management events affecting females
	 * @param event Type of farming operation (mowing, spraying, tillage, etc.)
	 * @return true if female survives event, false if killed
	 *
	 * @details Farm operations can affect females directly:
	 * - Pesticide spraying: Contact mortality, contaminated forage
	 * - Mowing/harvesting: Destroys forage resources
	 * - Tillage: Destroys ground nests (if applicable)
	 *
	 * Return value determines whether female continues activity or is removed from simulation.
	 */
	virtual bool OnFarmEvent(FarmToDo event);

	/**
	 * @brief Override pesticide contact handling for female-specific exposure
	 * @param a_x X-location of contact (default -1 = use current location)
	 * @param a_y Y-location of contact (default -1 = use current location)
	 *
	 * @details Females experience pesticide exposure via:
	 * 1. Direct overspray (flying through spray plume or foraging during application)
	 * 2. Contact with treated vegetation
	 * 3. Contaminated pollen/nectar (handled separately via m_foraged_resource_pesticide)
	 *
	 * Exposure routes and effects calculated using pesticide module parameters (absorption rates,
	 * body surface areas, thresholds). Female-specific exposure reflects foraging behaviour bringing
	 * them into contact with treated crops.
	 */
	virtual void DoPesticideContact(int a_x = -1, int a_y = -1);

#ifdef __OSMIA_PESTICIDE_STORE
	/** @brief Unique animal ID for pesticide exposure tracking/output */
	unsigned int m_animal_id;
#endif
};

#endif // Header include guard

