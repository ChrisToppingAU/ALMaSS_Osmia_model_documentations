/**
 * @class Osmia_Egg
 * @brief First life stage - egg developing within sealed nest cell
 * @extends Osmia_Base
 * 
 * @details Osmia_Egg represents the egg stage from laying through hatching. Development is driven
 * by degree-day accumulation above the threshold temperature. Eggs experience daily mortality risk
 * and may be parasitised if the cell was compromised during the vulnerable period whilst open.
 * 
 * @par Biological Foundation
 * *O. bicornis* eggs are laid on top of provision masses in sealed cells. Development time depends
 * strongly on temperature - warmer conditions accelerate hatching, cooler temperatures slow or halt
 * development. Laboratory studies (Giejdasz & Wilkaniec 2002, Radmacher & Strohm 2011) provide
 * baseline parameters, though implementation values are calibrated for field realism.
 * 
 * @par Development Model
 * Uses simple degree-day accumulation: each day above threshold temperature adds (T - T_threshold)
 * to m_AgeDegrees. When m_AgeDegrees reaches m_OsmiaEggDevelTotalDD (86 DD), egg hatches and
 * individual transitions to Osmia_Larva.
 * 
 * @par Mortality
 * Daily mortality probability applied each day (m_DailyDevelopmentMortEggs = 0.0014). Represents
 * combined effects of desiccation, fungal infection, temperature extremes, and developmental failures.
 * Mortality is temperature-independent in implementation despite some evidence of temperature effects.
 * 
 * @par Sex Determination
 * Sex is determined at laying (set by parent female based on sex allocation decisions). Female eggs
 * are fertilised (diploid), male eggs unfertilised (haploid). Sex persists through all life stages
 * though males are not explicitly modelled as adults.
 * 
 * @par Difference from Formal Model
 * The formal model specified development parameters; implementation adjusted them substantially for
 * field realism (see Osmia_Base parameter documentation). Mortality implementation matches formal
 * model exactly.
 * 
 * @see Osmia_Base for development parameters
 * @see st_Develop() for degree-day accumulation logic
 * @see st_Hatch() for transition to larval stage
 */
class Osmia_Egg : public Osmia_Base
{
protected:
	/**
	 * @var m_AgeDegrees
	 * @brief Accumulated degree-days toward hatching threshold
	 * @details Incremented daily by (T_today - T_threshold) when temperature exceeds threshold.
	 * When m_AgeDegrees >= m_OsmiaEggDevelTotalDD, egg hatches.
	 * 
	 * @par Implementation Note
	 * Initialized to 0.0 at egg creation. Can be queried for monitoring development progress.
	 * Persists through metamorphosis - larvae inherit m_AgeDegrees to continue development tracking.
	 */
	double m_AgeDegrees = 0.0;
	
	/**
	 * @var m_Sex
	 * @brief Sex of individual (true = female = fertilised egg)
	 * @details Determined at laying by parent female's sex allocation algorithm. Female eggs are
	 * fertilised (diploid, sex=true), male eggs unfertilised (haploid, sex=false). Sex persists
	 * through all life stages.
	 * 
	 * @par Biological Basis
	 * Hymenoptera have haplodiploid sex determination: fertilised eggs develop as females, unfertilised
	 * as males. Females control fertilisation during egg laying, enabling strategic sex allocation.
	 * *O. bicornis* typically places female eggs (larger provisions) deep in nest, male eggs (smaller
	 * provisions) near entrance.
	 */
	bool m_Sex;
	
	/**
	 * @var m_StageAge
	 * @brief Age in days when current life stage was entered
	 * @details Records absolute age (m_Age) at stage transition. Used to calculate stage duration
	 * for output and validation. For eggs, m_StageAge equals 0 (stage entered at birth).
	 */
	int m_StageAge;
	
	/**
	 * @var m_egg_pest_mortality
	 * @brief Cumulative pesticide mortality probability for this egg
	 * @details Tracks pesticide exposure effects when pesticide module is active. Represents both
	 * maternal transfer of pesticides via contaminated provisions and direct exposure if nest is
	 * sprayed. Used primarily for pesticide scenario testing.
	 */
	double m_egg_pest_mortality;

public:
	/**
	 * @brief Constructor for new egg object
	 * @param data Initialization data including parent info, nest location, sex, provision mass
	 * @details Initializes egg attributes, sets m_AgeDegrees to 0, records sex, links to nest.
	 * Egg begins development immediately after creation.
	 */
	Osmia_Egg(struct_Osmia* data);
	
	/**
	 * @brief Reinitialize egg object from object pool
	 * @param data New initialization data
	 * @details Resets all egg-specific attributes for object reuse, avoiding allocation overhead.
	 */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_Egg();
	
	/**
	 * @brief Main step function executing egg behaviour
	 * @details Calls st_Develop to accumulate degree-days and test for hatching. Applies daily
	 * mortality test. Transitions to Osmia_Larva upon hatching or to death upon mortality.
	 */
	virtual void Step(void);
	
	/** @brief Get accumulated degree-days */
	double GetAgeDegrees() { return m_AgeDegrees; }
	
	/** @brief Set accumulated degree-days (used during object reinitialization) */
	void SetAgeDegrees(unsigned a_agedegrees) { m_AgeDegrees = a_agedegrees; }

protected:
	/**
	 * @brief Development state - accumulate degree-days toward hatching
	 * @return Next state (st_Hatch if threshold reached, st_Develop to continue, or st_Die)
	 * 
	 * @details Each day:
	 * 1. Check if temperature exceeds threshold (m_OsmiaEggDevelThreshold)
	 * 2. If yes: add (T - threshold) to m_AgeDegrees
	 * 3. Check if m_AgeDegrees >= m_OsmiaEggDevelTotalDD
	 * 4. If yes: return st_Hatch
	 * 5. Apply daily mortality test
	 * 6. If survives: return st_Develop (continue development)
	 * 
	 * @par Temperature Handling
	 * Days below threshold add zero DD - development is suspended but mortality still applies.
	 * This can extend absolute calendar time to hatching if prolonged cool periods occur.
	 * 
	 * @par Difference from Formal Model
	 * **IMPLEMENTATION MATCH** - Logic follows formal model specification exactly. Parameter values
	 * differ (see Osmia_Base documentation) but algorithm structure is identical.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Transition state - metamorphose from egg to larva
	 * @return Next state (always returns st_Die as this object is deleted)
	 * 
	 * @details Hatching sequence:
	 * 1. Signal population manager to create new Osmia_Larva object
	 * 2. Transfer essential attributes (sex, mass, nest pointer, age, parasitism status, m_AgeDegrees)
	 * 3. Population manager updates nest's cell pointer to new larva object
	 * 4. This egg object is deleted
	 * 5. Larva begins its developmental cycle
	 * 
	 * @par Implementation Pattern
	 * This metamorphosis pattern (delete old object, create new object of next class, transfer attributes)
	 * is used at all stage transitions. It maintains clean separation between life stages whilst
	 * preserving individual identity through attribute transfer.
	 * 
	 * @par Biological Accuracy
	 * Represents the ecological concept of "hatching" - egg case breaks, first instar larva emerges
	 * and begins feeding on provision mass. In model terms, this is object type change without loss
	 * of individual continuity.
	 */
	virtual TTypeOfOsmiaState st_Hatch(void);
	
	/**
	 * @brief Daily mortality test for eggs
	 * @return true if dies, false if survives
	 * @details Simple probabilistic test: generates uniform random number [0,1], compares to
	 * m_DailyDevelopmentMortEggs (0.0014). Independent of temperature, age, or provision mass.
	 * 
	 * @par Biological Interpretation
	 * Represents average daily hazard from all mortality sources. Temperature-independence is
	 * simplification due to insufficient data on temperature-mortality relationships. Conservative
	 * choice given uncertainty.
	 */
	virtual bool DailyMortality() { 
		if (g_rand_uni_fnc() < m_DailyDevelopmentMortEggs) return true; 
		else return false; 
	}
};

/**
 * @class Osmia_Larva
 * @brief Feeding larval stage consuming provision mass
 * @extends Osmia_Egg
 * 
 * @details Osmia_Larva represents the actively feeding larval stage from hatching through cocoon
 * spinning. Larvae consume the provision mass left by their mother, growing through multiple instars
 * before entering the prepupal stage. Development continues via degree-day accumulation with
 * temperature threshold and requirements distinct from egg stage.
 * 
 * @par Biological Foundation
 * *O. bicornis* larvae progress through approximately 4-5 instars, consuming the entire provision
 * mass over 3-4 weeks (temperature dependent). Larvae are relatively sedentary within their sealed
 * cells, protected from most environmental hazards. Provision quality (pollen source diversity,
 * nutrient content) affects larval growth and survival, though this is not explicitly modelled beyond
 * mass effects.
 * 
 * @par Development Model
 * Continues degree-day accumulation from egg stage: m_AgeDegrees carries forward and increments by
 * (T - T_larva_threshold) each day. When m_AgeDegrees reaches m_OsmiaEggDevelTotalDD +
 * m_OsmiaLarvaDevelTotalDD, larva spins cocoon and transitions to prepupal stage.
 * 
 * @par Mass and Provisioning
 * Larval mass (m_Mass) represents the provision mass available in cell, set when egg was provisioned.
 * This determines final adult size via conversion equation. Larvae don't actively "eat" provision
 * (no depletion simulation) - mass is simply carried through as determinant of adult size.
 * 
 * @par Difference from Formal Model
 * Larval development threshold calibrated from 8.5°C to 4.5°C (see Osmia_Base). Mortality matches
 * formal model exactly. The simplified mass handling (no explicit feeding simulation) is consistent
 * with formal model's focus on outcome (adult mass) rather than process (feeding dynamics).
 * 
 * @see Osmia_Base for larval development parameters
 * @see st_Develop() for development logic
 * @see st_Prepupate() for transition to prepupal stage
 */
class Osmia_Larva : public Osmia_Egg
{
public:
	/**
	 * @brief Constructor for new larva (created during egg hatching)
	 * @param data Initialization data transferred from egg plus larval-specific attributes
	 * @details Inherits m_AgeDegrees from egg, increments m_StageAge to record larval entry,
	 * maintains sex, mass, and nest linkage. Ready to continue degree-day accumulation toward
	 * prepupal threshold.
	 */
	Osmia_Larva(struct_Osmia* data);
	
	/**
	 * @brief Reinitialize larva object from object pool
	 * @param data New initialization data
	 */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_Larva();
	
	/**
	 * @brief Main step function executing larval behaviour
	 * @details Very similar to egg step: calls st_Develop for degree-day accumulation, tests daily
	 * mortality, transitions to prepupa when threshold reached.
	 */
	virtual void Step(void);

protected:
	/**
	 * @brief Development state - accumulate degree-days toward prepupation
	 * @return Next state (st_Prepupate if threshold reached, st_Develop to continue, or st_Die)
	 * 
	 * @details Development logic parallel to egg stage but with larval parameters:
	 * 1. Add DD if T > m_OsmiaLarvaDevelThreshold (4.5°C)
	 * 2. Check if total DD >= (egg DD + larva DD requirements)
	 * 3. Transition to prepupa or continue development
	 * 4. Apply daily mortality test
	 * 
	 * @par Biological Process
	 * During larval development, individual progresses through instars (not explicitly modelled),
	 * consuming provision mass and growing. Final instar larva spins silk cocoon within cell,
	 * marking transition to prepupal stage.
	 * 
	 * @par Implementation Note
	 * The cumulative DD check (egg + larva requirements) means early hatching (warm egg conditions)
	 * can lead to slightly shorter larval duration and vice versa, creating realistic thermal
	 * integration across stages.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Transition state - metamorphose from larva to prepupa
	 * @return Next state (st_Die as larva object is deleted)
	 * 
	 * @details Prepupation sequence:
	 * 1. Signal population manager to create Osmia_Prepupa object
	 * 2. Transfer attributes (sex, mass, nest, age, DD, parasitism)
	 * 3. Update nest cell pointer to new prepupa
	 * 4. Delete larva object
	 * 5. Prepupa begins time-based development (not DD-based)
	 * 
	 * @par Biological Transition
	 * Final instar larva completes feeding, voids gut, spins cocoon, and enters prepupal diapause.
	 * This is a major physiological transition from feeding to dormancy, preparing for metamorphosis
	 * to adult form.
	 */
	virtual TTypeOfOsmiaState st_Prepupate(void);
	
	/**
	 * @brief Daily mortality test for larvae
	 * @return true if dies, false if survives
	 * @details Uses m_DailyDevelopmentMortLarvae (0.0014), identical to egg mortality. Temperature-
	 * independent despite some evidence that provision quality (temperature-affected) influences
	 * survival.
	 */
	virtual bool DailyMortality() { 
		if (g_rand_uni_fnc() < m_DailyDevelopmentMortLarvae) return true; 
		else return false; 
	}
};

/**
 * @class Osmia_Prepupa
 * @brief Prepupal diapause stage in cocoon
 * @extends Osmia_Larva
 * 
 * @details Osmia_Prepupa represents the prepupal dormancy period following cocoon spinning. This
 * stage uses time-based rather than degree-day-based development, with non-linear temperature
 * relationships. Prepupae are relatively invulnerable, having lowest mortality of any stage.
 * 
 * @par Biological Foundation
 * The prepupal stage is characterized by arrested development (diapause) lasting 1-3 months depending
 * on temperature and photoperiod cues. *O. bicornis* prepupae show optimal development at intermediate
 * temperatures (~22°C) with both lower and higher temperatures extending development time. This
 * non-linear response distinguishes prepupae from other stages' monotonic temperature relationships.
 * 
 * @par Development Model
 * **MAJOR DIFFERENCE FROM FORMAL MODEL**: Formal model specified quadratic temperature-development
 * function with 24.3-day optimum at 22°C. Implementation uses simpler time-based approach: base
 * duration 45 days with individual variation (±10%), plus temperature threshold effects rather than
 * continuous function. This pragmatic simplification reflects limited data for robust parameterization
 * of non-linear relationship.
 * 
 * @par Rationale for Simplified Model
 * Quadratic functions require precise parameterization of both optimal temperature and curvature.
 * Available data (Radmacher & Strohm 2011, Giejdasz & Fliszkiewicz 2016) show scatter making curve
 * fitting uncertain. Time-based approach with thresholds provides more stable model behaviour whilst
 * capturing key biology: prepupae take several weeks and respond to temperature extremes.
 * 
 * @par Mortality
 * Lowest of all stages (0.003 daily probability). Well-supported empirically - cocooned prepupae
 * are protected and physiologically inactive, minimizing vulnerability.
 */
class Osmia_Prepupa : public Osmia_Larva
{
public:
	/** @brief Constructor transferring from larva to prepupa */
	Osmia_Prepupa(struct_Osmia* data);
	
	/** @brief Reinitialize prepupa from object pool */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_Prepupa();
	
	/** @brief Main step function */
	virtual void Step(void);

protected:
	/**
	 * @brief Development state - time-based progression toward pupation
	 * @return Next state (st_Pupate when duration complete, st_Develop to continue, or st_Die)
	 * 
	 * @details Different development logic from egg/larva:
	 * 1. Check if prepupal age >= m_myOsmiaPrepupaDevelTotalDays
	 * 2. Temperature affects this duration via threshold checks (not continuous function)
	 * 3. Apply daily mortality test (very low)
	 * 4. Transition to pupa when time threshold reached
	 * 
	 * @par Temperature Effects
	 * Rather than quadratic function, implementation uses temperature thresholds:
	 * - Above prewintering threshold (15°C): development proceeds
	 * - Below threshold: development suspended
	 * This creates appropriate seasonal timing without complex non-linear function.
	 * 
	 * @par Individual Variation
	 * Each prepupa gets individual m_myOsmiaPrepupaDevelTotalDays drawn from uniform distribution
	 * (base ± 10%), creating realistic spread in prepupal durations even under identical temperatures.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Transition state - metamorphose from prepupa to pupa
	 * @return Next state (st_Die as prepupa deleted)
	 * @details Standard metamorphosis: create Osmia_Pupa, transfer attributes, update nest, delete prepupa.
	 * Biologically represents transition from dormant prepupa to active metamorphic pupa.
	 */
	virtual TTypeOfOsmiaState st_Pupate(void);
	
	/**
	 * @brief Daily mortality test for prepupae
	 * @return true if dies, false if survives
	 * @details Uses m_DailyDevelopmentMortPrepupae (0.003) - higher than egg/larva but still very low.
	 * Cocooned prepupae are well-protected.
	 */
	virtual bool DailyMortality() { 
		if (g_rand_uni_fnc() < m_DailyDevelopmentMortPrepupae) return true; 
		else return false; 
	}
	
	/**
	 * @var m_myOsmiaPrepupaDevelTotalDays
	 * @brief Individual-specific prepupal development duration (days)
	 * @details Drawn at prepupa creation from uniform distribution: base value (45 days) ±10%.
	 * Creates individual variation in development timing. Value remains constant for this individual's
	 * prepupal stage.
	 * 
	 * @par Biological Basis
	 * Prepupal duration varies substantially between individuals even under identical conditions,
	 * reflecting genetic variation, maternal effects, and provision quality effects not explicitly
	 * captured otherwise. This stochastic variation prevents unrealistic synchronous pupation.
	 */
	double m_myOsmiaPrepupaDevelTotalDays;
};

/**
 * @class Osmia_Pupa
 * @brief Pupal metamorphosis stage
 * @extends Osmia_Prepupa
 * 
 * @details Osmia_Pupa represents active metamorphosis from larval to adult form. Returns to degree-day
 * development model (like egg/larva) after prepupal time-based development. Pupal stage ends with
 * eclosion - adult emerges from pupal exuviae within cocoon, entering overwintering adult stage.
 * 
 * @par Biological Foundation
 * During pupation, tissues reorganize from larval to adult configuration through histolysis and
 * histogenesis. This metabolically expensive process is temperature-sensitive - higher temperatures
 * accelerate metamorphosis. Pupae remain within cocoon, protected but vulnerable to extreme temperatures.
 * 
 * @par Development Model
 * Returns to degree-day accumulation after prepupal time-based development. However, m_AgeDegrees
 * is reset at pupation (not carried forward from earlier stages) because pupal DD requirements are
 * independent of earlier development history. Accumulates (T - threshold) daily until reaching
 * m_OsmiaPupaDevelTotalDD, then transitions to overwintering adult.
 * 
 * @par Calibration
 * **MAJOR DIFFERENCE**: Pupal parameters show largest calibration adjustment in model. Formal model:
 * 272.3 DD with 13.2°C threshold. Implementation: 570 DD with 1.1°C threshold. This compensatory
 * adjustment maintains realistic development timing under field temperatures whilst preventing
 * developmental failure. Code comment explicitly notes: "changed from 13.2 to prevent pupal death".
 * 
 * @par Mortality
 * Very low (0.003 daily), same as prepupa. Well-protected within cocoon. Most pupal mortality
 * reflects earlier problems (insufficient provisioning) manifesting during energy-intensive metamorphosis.
 * 
 * @see Osmia_Base for detailed discussion of pupal parameter calibration
 */
class Osmia_Pupa : public Osmia_Prepupa
{
public:
	/** @brief Constructor transferring from prepupa */
	Osmia_Pupa(struct_Osmia* data);
	
	/** @brief Reinitialize pupa from object pool */
	virtual void ReInit(struct_Osmia* data);
	
	/** @brief Destructor */
	virtual ~Osmia_Pupa();
	
	/** @brief Main step function */
	virtual void Step(void);

protected:
	/**
	 * @brief Development state - accumulate degree-days toward eclosion
	 * @return Next state (st_Emerge when DD threshold reached, st_Develop to continue, or st_Die)
	 * 
	 * @details Standard degree-day logic:
	 * 1. If T > m_OsmiaPupaDevelThreshold (1.1°C): add (T - 1.1) to m_AgeDegrees
	 * 2. When m_AgeDegrees >= m_OsmiaPupaDevelTotalDD (570 DD): transition to adult
	 * 3. Apply daily mortality test
	 * 
	 * @par Temperature Response
	 * Very low threshold (1.1°C) means pupal development proceeds across full summer temperature
	 * range. High total DD requirement (570) compensates, maintaining appropriate absolute duration.
	 * This parameterization essential for preventing developmental failures noted during calibration.
	 * 
	 * @par Biological Process
	 * Pupal stage is active metamorphosis. Adult structures (wings, legs, reproductive organs) form
	 * from imaginal discs whilst larval tissues are remodelled. Process requires substantial energy
	 * and time, reflected in high DD requirement.
	 */
	virtual TTypeOfOsmiaState st_Develop(void);
	
	/**
	 * @brief Transition state - eclosion from pupa to adult-in-cocoon
	 * @return Next state (st_Die as pupa deleted)
	 * 
	 * @details Eclosion sequence:
	 * 1. Create Osmia_InCocoon object (overwintering adult)
	 * 2. Transfer attributes including final adult mass (from provision mass conversion)
	 * 3. Update nest pointer
	 * 4. Delete pupa object
	 * 5. Adult-in-cocoon begins overwintering phase
	 * 
	 * @par Biological Transition
	 * Adult bee emerges from pupal cuticle but remains within protective cocoon. Will not leave nest
	 * until spring emergence conditions are met (diapause completion + warming). This is key
	 * overwintering strategy - adults overwinter in cocoons rather than as pupae.
	 */
	virtual TTypeOfOsmiaState st_Emerge(void);
	
	/**
	 * @brief Daily mortality test for pupae
	 * @return true if dies, false if survives
	 * @details Uses m_DailyDevelopmentMortPupae (0.003). Low mortality reflects protection of cocooned
	 * pupae. Most failures occur if larva was inadequately provisioned, causing pupal metamorphosis
	 * to fail (not explicitly modelled - assumed rare).
	 */
	virtual bool DailyMortality() { 
		if (g_rand_uni_fnc() < m_DailyDevelopmentMortPupae) return true; 
		else return false; 
	}
};

#endif // Include guard for header file continuation
