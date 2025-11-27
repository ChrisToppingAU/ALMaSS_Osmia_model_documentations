/*
*******************************************************************************************************
Copyright (c) 2019, Christopher John Topping, Aarhus University
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the
following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************************************************
*/

/**
 * @file Osmia.h
 * @brief Agent-based implementation of *Osmia bicornis* (red mason bee) life cycle
 * 
 * @details This header file defines the complete class hierarchy for simulating *Osmia bicornis*
 * populations within the ALMaSS framework. The implementation follows the formal model described in
 * Ziółkowska et al. (2025), covering all life stages from egg through adult, including temperature-driven
 * development, overwintering physiology, foraging behaviour, and nest provisioning.
 * 
 * @par Biological Foundation
 * The model is parameterised primarily from laboratory studies by Radmacher and Strohm (2011), 
 * Giejdasz and Wilkaniec (2002), and Giejdasz and Fliszkiewicz (2016) for developmental rates and 
 * thresholds. Foraging behaviour draws from field observations by Seidelmann (2006), whilst 
 * overwintering mortality relationships come from Sgolastra et al. (2011) working with 
 * *O. lignaria*.
 * 
 * @par Implementation Approach
 * The model implements a stage-structured agent-based approach where each individual progresses through
 * discrete life stages (egg, larva, prepupa, pupa, overwintering adult, active adult). Development is
 * primarily temperature-driven using degree-day accumulation, with mortality applied as daily 
 * probabilities at each stage. Spatial behaviour emerges from individual movement and foraging 
 * decisions based on local resource availability and distance constraints.
 * 
 * @par Key Design Decisions
 * - Males are not explicitly modelled; reproductive success focuses on female provisioning and egg production
 * - Prepupal development uses time-based rather than degree-day approach due to non-linear temperature response
 * - Foraging employs a detailed spatial mask for efficient resource searches without repeated distance calculations
 * - Nests are modelled as linear structures with sequentially provisioned cells
 * 
 * @author Original implementation: Christopher J. Topping
 * @author Enhanced documentation: [Author Name]
 * @date Original: August 2019
 * @date Enhanced: 2025
 * @ingroup Osmia_Model
 * 
 * @see Osmia.cpp for implementation details
 * @see Osmia_Population_Manager.h for population-level management
 * @see Ziółkowska et al. (2025) Food and Ecological Systems Modelling Journal for formal model specification
 */

//---------------------------------------------------------------------------
#ifndef OsmiaH
#define OsmiaH
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include <forward_list>

class Osmia_Population_Manager;
class OsmiaParasitoid_Population_Manager;
class Osmia_Nest_Manager;
class PollenMap_centroidbased;
class probability_distribution;
class Osmia_Egg;
class Osmia_Female;
class struct_Osmia;
class Osmia_Base;

//------------------------------------------------------------------------------
/**
 * Used for the population manager's list of Osmia
 */
//typedef vector<Osmia*> TListOfOsmia;
//---------------------------------------------------------------------------

/**
 * @def __OSMIA_DIST_SIZE
 * @brief Size of pre-calculated distribution arrays for movement probabilities
 * @details This size determines the resolution of probability distributions used for dispersal
 * and foraging movements. Larger values provide finer resolution but increase memory usage.
 * Value of 10,000 provides sufficient precision for beta distributions with shape parameters
 * typically in range 1-10.
 */
#define __OSMIA_DIST_SIZE 10000

/**
 * @var cfg_OsmiaTotalCocoonMassLoss
 * @brief Total mass loss from first to last cocoon per nest (female cocoon mass units)
 * @details Default: 15.0 mg
 * 
 * @par Empirical Basis
 * Based on observations from Seidelmann (2010) showing that female *O. bicornis* exhibit declining
 * cocoon masses from first to last offspring within a nest, reflecting progressive depletion of 
 * maternal resources. The value represents average total decline across a complete nest.
 * 
 * @par Biological Interpretation
 * This progressive mass loss reflects the declining foraging efficiency and accumulated physiological
 * costs as females age. Later offspring receive slightly less provisioning, potentially affecting
 * their survival and future reproductive success.
 * 
 * @par Implementation Note
 * The mass loss is distributed across nest cells with added stochastic variation (see
 * cfg_OsmiaTotalCocoonMassLossRange). This creates realistic within-nest variation in offspring
 * condition.
 * 
 * @par Uncertainty
 * MEDIUM - Individual variation in foraging success and longevity creates substantial variation
 * around this mean value. Field conditions may produce different patterns than semi-controlled
 * observation studies.
 */
static CfgFloat cfg_OsmiaTotalCocoonMassLoss("OSMIATOTALCOCOONMASSLOSS", CFG_CUSTOM, 15.0);

/**
 * @var cfg_OsmiaTotalCocoonMassLossRange
 * @brief Stochastic range around total cocoon mass loss between nests (female cocoon mass units)
 * @details Default: 5.0 mg
 * 
 * @par Biological Rationale
 * Individual females vary in foraging ability, nest location quality, and accumulated wear, 
 * creating variation in the magnitude of progressive mass loss to offspring. This parameter
 * captures between-female variation in resource provisioning patterns.
 * 
 * @par Implementation
 * Applied as ±range around the mean total mass loss, creating a uniform distribution of 
 * possible mass loss trajectories across females in the population.
 * 
 * @par Uncertainty
 * LOW - The existence of substantial variation is well established; this value provides
 * reasonable spread without extreme outliers.
 */
static CfgFloat cfg_OsmiaTotalCocoonMassLossRange("OSMIATOTALCOCOONMASSLOSSRANGE", CFG_CUSTOM, 5.0);

/**
 * @enum TTypeOfOsmiaState
 * @brief Behavioural states governing *Osmia bicornis* agent decisions
 * 
 * @details This enumeration defines the discrete behavioural states that structure the decision-making
 * of *Osmia* agents. Each state represents a distinct mode of behaviour with specific rules and 
 * possible transitions. The state machine approach provides clear separation of behavioural logic
 * and facilitates debugging and model extension.
 * 
 * @par Implementation Pattern
 * States are implemented via virtual functions (st_StateName) that return the next state to transition
 * to. This allows state-specific behaviour to be defined in derived classes whilst maintaining
 * a common control flow structure.
 * 
 * @par State Transition Logic
 * Most states follow the pattern: perform behaviour → check conditions → return next state or toOsmias_Die.
 * The population manager calls Step() repeatedly until all agents return a terminal state.
 */
enum TTypeOfOsmiaState
{
	/** @brief Initial state upon object creation; performs setup and transitions to first active state */
	toOsmias_InitialState = 0,
	
	/** @brief Active development state; accumulates degree-days or time towards stage transition */
	toOsmias_Develop,
	
	/** @brief Transition state for metamorphosis to next life stage; handles object type conversion */
	toOsmias_NextStage,
	
	/** @brief Dispersal state for adult females seeking nesting locations beyond their natal area */
	toOsmias_Disperse,
	
	/** @brief Active provisioning state; foraging for pollen and nectar to stock nest cells */
	toOsmias_NestProvisioning,
	
	/** @brief Reproductive decision-making state; includes nest finding, sex allocation, and egg laying */
	toOsmias_ReproductiveBehaviour,
	
	/** @brief Post-emergence state before initiating reproduction; includes maturation and mating */
	toOsmias_Emerged,
	
	/** @brief Terminal state; agent removed from simulation */
	toOsmias_Die
};

/**
 * @enum TTypeOfOsmiaParasitoids
 * @brief Classification of parasitoid types affecting *Osmia bicornis* eggs and larvae
 * 
 * @details This enumeration defines the possible parasitism outcomes for developing *Osmia* individuals
 * within nests. Different parasitoid types have distinct attack probabilities and timing, affecting
 * host survival differently.
 * 
 * @par Biological Background
 * *O. bicornis* nests are susceptible to various natural enemies including: (1) bombyliid flies that
 * enter open nest cells and lay eggs on or near the provision mass, (2) cleptoparasitic bees that
 * steal provisions, and (3) other parasitoids. Parasitism risk increases with nest cell open time.
 * 
 * @par Implementation
 * The enumeration uses unsigned type because values may be used as array indices for parasitoid-specific
 * parameters. Each parasitised individual is marked with one parasitoid type; multiple parasitism is
 * not currently modelled.
 * 
 * @par Difference from Formal Model
 * The formal model described conceptual parasitism mechanisms but left specific parameterisation to
 * the implementation. This enum structure supports future addition of detailed parasitoid dynamics.
 * 
 * @see CalcParaistised() for parasitism determination logic
 */
enum class TTypeOfOsmiaParasitoids : unsigned
{
	/** @brief Egg/larva develops normally without parasitism */
	topara_Unparasitised = 0,
	
	/** @brief Parasitised by bombyliid fly; typically lethal to host */
	topara_Bombylid,
	
	/** @brief Provisions stolen by cleptoparasitic bee; host starves */
	topara_Cleptoparasite,
	
	/** @brief Placeholder for future parasitoid types */
	topara_foobar
};

/**
 * @class OsmiaForageMask
 * @brief Pre-calculated spatial search mask for efficient resource location
 * 
 * @details This class provides a pre-computed mask of spatial offsets that can be iterated through
 * without repeated distance and direction calculations. The mask defines a series of concentric 
 * distance bands with eight cardinal/intercardinal directions at each distance.
 * 
 * @par Performance Rationale
 * Resource searches are called frequently (multiple times per foraging bee per day), making 
 * computational efficiency critical. Pre-calculating offset patterns eliminates trigonometric
 * operations within search loops, substantially reducing execution time in landscape-scale simulations.
 * 
 * @par Biological Interpretation
 * The radial search pattern reflects observations that *Osmia* females tend to forage progressively
 * further from their nest if nearby resources are depleted or of poor quality. The 8-direction
 * structure approximates the actual multi-directional search behaviour whilst keeping memory usage
 * tractable.
 * 
 * @par Implementation Details
 * The mask array holds [distance_step][direction][x or y offset] values. Step size can be configured
 * to match landscape resolution. Typical usage iterates through distances (nearest first) and 
 * directions, testing each location for resource availability.
 * 
 * @see OsmiaForageMaskDetailed for higher-resolution variant
 */
class OsmiaForageMask
{
public:
	/** 
	 * @brief Three-dimensional array: [20 distances][8 directions][2 coordinates (x,y)]
	 * @details Structure stores integer offsets from a centre point for efficient spatial searches.
	 * Twenty distance steps provide coverage to typical foraging ranges; eight directions balance
	 * coverage with memory usage.
	 */
	int m_mask[20][8][2];
	
	/** 
	 * @brief Step size in landscape units between successive distance rings
	 * @details Determines the granularity of the search. Smaller steps provide finer coverage but
	 * require more iterations. Typical value is 1-2 landscape grid cells.
	 */
	int m_step;
	
	/** 
	 * @brief Squared step size for distance calculations
	 * @details Pre-computed to avoid repeated multiplication in distance comparisons. Used when
	 * checking if a location falls within the current search radius.
	 */
	int m_step2;
	
	/**
	 * @brief Constructor initialising the spatial offset mask
	 * @details Calculates and stores offset values for all distance-direction combinations.
	 * Called once during population manager initialisation to avoid repeated computation.
	 */
	OsmiaForageMask();
};

/**
 * @class OsmiaForageMaskDetailed
 * @brief High-resolution spatial search mask for detailed resource assessment
 * 
 * @details This variant provides finer-grained spatial coverage than OsmiaForageMask, storing all
 * offsets as a sequential vector rather than distance-direction arrays. Used when detailed spatial
 * analysis is required, trading increased memory usage for improved coverage.
 * 
 * @par Usage Context
 * Employed primarily for pollen resource assessment where identifying all cells within foraging range
 * is more important than the progressive search strategy. The vector structure simplifies iteration
 * when the search order is less critical.
 * 
 * @par Difference from OsmiaForageMask
 * Whilst OsmiaForageMask provides coarse distance bands with 8 directions per band, this class stores
 * all locations within maximum distance as a flat vector. This supports different search algorithms
 * (e.g., parallel assessment of all available resources vs. incremental outward search).
 */
class OsmiaForageMaskDetailed
{
public:
	/** 
	 * @brief Vector storing all spatial offsets within maximum distance
	 * @details Each APoint contains x,y coordinates relative to the search centre. Points are stored
	 * in an order facilitating efficient iteration, typically organised by distance from centre.
	 */
	vector<APoint> m_mask;
	
	/** 
	 * @brief Step size in landscape units between sampled locations
	 * @details Controls search resolution. Step of 1 samples every grid cell; larger steps reduce
	 * computational load but may miss small resource patches.
	 */
	int m_step;
	
	/** 
	 * @brief Maximum search distance in landscape units
	 * @details Defines the outer boundary of the search mask. Typically set to match species-specific
	 * foraging range constraints from homing distance data.
	 */
	int m_maxdistance;
	
	/**
	 * @brief Constructor creating detailed mask with specified resolution and range
	 * @param a_step Step size between sampled points
	 * @param a_maxdistance Maximum radius of search mask
	 * @details Generates all offset points within the specified distance, storing them in the mask vector.
	 */
	OsmiaForageMaskDetailed(int a_step, int a_maxdistance);
};

/**
 * @class OsmiaNestData
 * @brief Data structure recording nest contents and provisioning status
 * 
 * @details Simple container class tracking the current state of a nest under construction. Used
 * primarily for testing and validation to compare intended provisioning plans with actual outcomes.
 * 
 * @par Usage
 * When compiled with __OSMIATESTING defined, instances record each female's target vs. achieved
 * provisioning. This enables post-simulation analysis of how environmental constraints (resource
 * availability, weather, mortality) affect reproductive success.
 */
class OsmiaNestData
{
public:
	/** @brief Number of eggs currently in the nest */
	int m_no_eggs;
	
	/** @brief Number of female eggs in the nest (for sex ratio analysis) */
	int m_no_females;
	
	/** @brief Vector recording provision mass (mg) for each nest cell in sequential order */
	vector<double> m_cell_provision;
};

/**
 * @class Osmia_Nest
 * @brief Container representing a linear nest structure with sequentially provisioned cells
 * @extends TAnimal
 * 
 * @details The Osmia_Nest class models the physical nest as a linear sequence of brood cells, each
 * containing a single *Osmia* egg or developing larva. The class primarily serves as a container,
 * maintaining pointers to its contained individuals and providing thread-safe access via locks.
 * 
 * @par Biological Basis
 * *O. bicornis* nests are naturally linear, created in pre-existing cavities (beetle borings, hollow
 * stems, trap nests). Females provision cells sequentially from the back of the cavity forward, 
 * placing an egg on each provision mass before sealing the cell and beginning the next. This structure
 * is faithfully represented in the model.
 * 
 * @par Inheritance from TAnimal
 * Extending TAnimal provides spatial location (x,y coordinates and polygon reference) and potential
 * access to the ALMaSS Step mechanism. However, nests are currently passive containers; their state
 * changes only through actions by the containing eggs/larvae or the provisioning female.
 * 
 * @par Thread Safety
 * The nest lock (m_cell_lock) prevents race conditions when multiple females might access nest data
 * simultaneously in parallelised simulations. All nest modifications must acquire the lock first.
 * 
 * @par Implementation Note
 * The cells are stored as a forward_list rather than vector because cells are only added (never 
 * removed or accessed by index), and forward_list provides efficient insertion with minimal memory
 * overhead.
 */
class Osmia_Nest : public TAnimal
{
protected:
	/** 
	 * @brief X-coordinate of nest location in landscape grid
	 * @details Position is set at nest creation and remains constant. Determines which landscape
	 * features (vegetation types, elevation) affect the nest environment.
	 */
	int m_x;
	
	/** 
	 * @brief Y-coordinate of nest location in landscape grid
	 * @details Used in conjunction with m_x for spatial queries and distance calculations.
	 */
	int m_y;
	
	/** 
	 * @brief Reference to landscape polygon containing the nest
	 * @details Links nest to specific landscape elements (habitat types, farm fields), enabling
	 * queries about local conditions and management events.
	 */
	int m_PolyRef;
	
	/** 
	 * @brief Forward-linked list of pointers to Osmia_Egg objects contained in nest cells
	 * @details Each element represents one provisioned cell with its egg or developing larva.
	 * Cells are added to the front of the list as the female provisions them sequentially.
	 * This maintains temporal order (newest cells at front).
	 */
	std::forward_list<TAnimal*>m_cells;
	
	/** 
	 * @brief OpenMP nested lock for thread-safe nest access
	 * @details Critical for preventing race conditions in parallel simulations where multiple agents
	 * might query or modify nest contents concurrently. Nested lock allows re-entrant access by the
	 * same thread if needed.
	 */
	omp_nest_lock_t* m_cell_lock;
	
	/** 
	 * @brief Static pointer to the single Osmia_Nest_Manager instance
	 * @details Provides all nests with access to population-level services (landscape queries,
	 * global configuration parameters, logging). Static ensures single shared manager.
	 */
	static Osmia_Nest_Manager* m_OurManager;
	
	/** 
	 * @brief Flag indicating whether nest is open for adding new cells
	 * @details Set to false when nest is sealed (female completes provisioning or dies). Prevents
	 * addition of new cells to abandoned nests. True whilst active female is provisioning.
	 */
	bool m_isOpen;
	
	/** 
	 * @brief Simulated micro-environmental variation in development timing (days)
	 * @details Represents aspect, exposure, and other micro-site effects causing individual nests to
	 * differ in thermal regime even at same location. Added as delay to emergence timing, creating
	 * realistic spread in emergence dates.
	 * 
	 * @par Biological Rationale
	 * Real nests experience thermal heterogeneity due to orientation (sunny vs. shaded), substrate
	 * type (wood vs. stems), and sheltering effects. This creates stochastic variation in development
	 * rates even for nearby nests, as documented in field emergence patterns showing 2-3 week spread.
	 * 
	 * @par Implementation
	 * Value assigned at nest creation from a distribution (typically normal or uniform). Applied as
	 * additive delay to thermal development calculations, simulating cooler micro-sites developing
	 * more slowly.
	 */
	int m_aspectdelay;

public:
	/**
	 * @brief Construct a new Osmia_Nest object at specified location
	 * @param a_x X-coordinate in landscape grid
	 * @param a_y Y-coordinate in landscape grid
	 * @param a_polyref Polygon reference for landscape context
	 * @param a_manager Pointer to nest population manager
	 * 
	 * @details Initialises nest at specified location, sets isOpen to true, calculates aspect delay,
	 * and creates the thread lock. The nest is ready to receive eggs from a provisioning female.
	 */
	Osmia_Nest(int a_x, int a_y, int a_polyref, Osmia_Nest_Manager* a_manager);
	
	/**
	 * @brief Destructor cleaning up nest resources
	 * @details Destroys the nest lock and frees associated memory. Contained egg/larva objects are
	 * managed separately by the population manager, so this destructor does not delete them.
	 */
	virtual ~Osmia_Nest(){
		omp_destroy_nest_lock(m_cell_lock);
		delete m_cell_lock;
	}
	
	/**
	 * @brief Acquire the nest lock for thread-safe access
	 * @details Call before any operation that reads or modifies m_cells. Blocks if another thread
	 * currently holds the lock. Nested lock allows same thread to re-acquire.
	 * 
	 * @par Usage Pattern
	 * Always pair with ReleaseCellLock() in same scope. Use RAII wrapper or ensure release even if
	 * exceptions occur. Typical pattern:
	 * @code
	 * SetCellLock();
	 * // ... critical section operations ...
	 * ReleaseCellLock();
	 * @endcode
	 */
	void SetCellLock(void) { omp_set_nest_lock(m_cell_lock); }
	
	/**
	 * @brief Release the nest lock after completing modifications
	 * @details Must be called after every SetCellLock() to prevent deadlocks. Allows waiting threads
	 * to proceed with their nest access.
	 */
	void ReleaseCellLock(void) { omp_unset_nest_lock(m_cell_lock); }
	
	/**
	 * @brief Add a cocoon to the nest (initialisation only)
	 * @param a_cocoon Pointer to overwintering individual
	 * 
	 * @details This method is used exclusively during simulation initialisation to populate nests with
	 * overwintering individuals from previous seasons. Not used during normal simulation runtime where
	 * eggs are added via AddEgg().
	 * 
	 * @par Implementation Note
	 * Uses push_front() to add to forward_list, which is O(1). Initialisation order (back to front or
	 * front to back) doesn't affect subsequent simulation as cocoons emerge based on temperature not
	 * position in nest.
	 */
	void AddCocoon (TAnimal* a_cocoon) {
		m_cells.push_front(a_cocoon);
	}
	
	/** 
	 * @brief Add a newly laid egg to the nest
	 * @param a_egg Pointer to Osmia_Egg object just created by provisioning female
	 * 
	 * @details Appends egg to the nest's cell list in a thread-safe manner. Called by Osmia_Female
	 * during egg laying after a cell has been fully provisioned.
	 * 
	 * @par Biological Timing
	 * Corresponds to the moment when the female seals a cell partition after placing an egg on the
	 * provision mass. At this point the egg begins development and becomes vulnerable to parasitism
	 * (if the cell seal is not perfect).
	 * 
	 * @par Thread Safety
	 * Method acquires cell lock internally, so calling code does not need to lock explicitly.
	 */
	void AddEgg (TAnimal* a_egg);
	
	/**
	 * @brief Replace a cell pointer during metamorphosis (egg→larva, larva→prepupa, etc.)
	 * @param a_old_ptr Pointer to object being replaced (about to be deleted)
	 * @param a_new_ptr Pointer to new life stage object
	 * 
	 * @details When an individual transitions between life stages, the old object is deleted and a new
	 * object of the appropriate class is created. This method updates the nest's cell list to point to
	 * the new object whilst maintaining cell order.
	 * 
	 * @par Implementation Detail
	 * Searches the forward_list for a_old_ptr and replaces it with a_new_ptr. Requires linear search
	 * through list, but list lengths are short (typically <15 cells) so performance is acceptable.
	 * 
	 * @par Thread Safety
	 * Caller must hold nest lock before calling this method to prevent concurrent modifications during
	 * pointer replacement.
	 */
	void ReplaceNestPointer(TAnimal* a_old_ptr, TAnimal* a_new_ptr);
	
	/**
	 * @brief Get count of cells currently in the nest
	 * @return Number of cells (eggs + developing larvae)
	 * 
	 * @details Returns the current size of m_cells list. Used for monitoring nest provisioning progress
	 * and for calculating parasitism risk (which increases with nest cell count).
	 * 
	 * @par Usage Note
	 * This count includes all cells added to date, including any that may have died. Dead cells are not
	 * actively removed from the list; they simply cease to Step and remain as inactive pointers until
	 * nest cleanup at season end.
	 */
	int GetCellCount();
	
	/**
	 * @brief Check if nest is open for adding new cells
	 * @return true if nest accepts new eggs, false if sealed or abandoned
	 * 
	 * @details Returns m_isOpen status. Closed nests reject new egg additions. Status is managed by
	 * the provisioning female, who sets it false when completing the nest or upon death.
	 */
	bool IsOpen() { return m_isOpen; }
	
	/**
	 * @brief Set nest open/closed status
	 * @param a_status true to open nest, false to close
	 * 
	 * @details Called by provisioning female when sealing the final nest cell or when abandoning a nest.
	 * Also may be set false by population manager during cleanup of nests belonging to dead females.
	 */
	void SetIsOpen(bool a_status) { m_isOpen = a_status; }
	
	/**
	 * @brief Get nest X-coordinate
	 * @return X position in landscape grid
	 */
	int GetX() { return m_x; }
	
	/**
	 * @brief Get nest Y-coordinate
	 * @return Y position in landscape grid
	 */
	int GetY() { return m_y; }
	
	/**
	 * @brief Get polygon reference for nest location
	 * @return Integer reference to landscape polygon
	 * 
	 * @details Used to query landscape manager for local conditions (habitat type, management events).
	 * Polygon reference remains constant for nest lifetime.
	 */
	int GetPolyRef() { return m_PolyRef; }
	
	/**
	 * @brief Get micro-environmental aspect delay
	 * @return Development delay in days
	 * 
	 * @details Returns the fixed delay assigned at nest creation, representing cooler/warmer micro-site
	 * effects. Used by overwintering individuals to adjust emergence timing.
	 */
	int GetAspectDelay() { return m_aspectdelay; }
};

/**
 * @class Osmia_Base
 * @brief Foundation class for all *Osmia bicornis* life stages
 * @extends TAnimal
 * 
 * @details Osmia_Base provides the common attributes, parameters, and methods shared across all
 * *O. bicornis* life stages from egg through adult female. The class holds static parameters for
 * development thresholds, mortality rates, and mass relationships that apply population-wide, plus
 * instance attributes tracking individual state (age, mass, current nest, parasitism status).
 * 
 * @par Biological Foundation
 * The class structure reflects the stage-structured life cycle of *O. bicornis* whilst maintaining
 * shared attributes that persist through metamorphosis (mass, sex, nest location, parasitism).
 * Static parameters ensure consistent application of population-level biology across all individuals.
 * 
 * @par Class Hierarchy Design
 * The inheritance chain (Osmia_Base → Osmia_Egg → Osmia_Larva → Osmia_Prepupa → Osmia_Pupa →
 * Osmia_InCocoon → Osmia_Female) allows progressive addition of stage-specific attributes and
 * behaviours whilst maintaining access to base functionality. During metamorphosis, objects are
 * deleted and recreated as the appropriate derived class, with key attributes copied forward.
 * 
 * @par Inheritance from TAnimal
 * Extends TAnimal to gain spatial location (m_Location_x, m_Location_y) and landscape access
 * (m_OurLandscape). This enables individuals to query local environmental conditions and respond
 * to landscape-level events (farming operations, weather).
 * 
 * @par Implementation Note - Static Parameters
 * Development and mortality parameters are static (shared across all instances) because they represent
 * population-level biology, not individual variation. This saves memory and ensures parameter consistency.
 * Individual variation emerges from stochastic processes (mortality tests, provision mass variation)
 * applied to these base parameters.
 * 
 * @par Difference from Formal Model
 * The formal model describes conceptual state variables and transitions. This implementation adds
 * necessary infrastructure (population manager pointers, thread-safe access, state machine
 * implementation) whilst faithfully representing the biological processes.
 */
class Osmia_Base : public TAnimal
{
protected:
	/**
	 * @var m_CurrentOState
	 * @brief Current behavioural state governing agent decisions
	 * @details State machine variable determining which behaviour method (st_Develop, st_Disperse,
	 * st_ReproductiveBehaviour, etc.) is active. Transitions between states occur based on
	 * developmental progress, environmental conditions, and mortality events.
	 */
	TTypeOfOsmiaState m_CurrentOState;
	
	/**
	 * @var m_Age
	 * @brief Chronological age in days since egg laying
	 * @details Incremented daily by the population manager. Used primarily for tracking life stage
	 * duration and for debugging/output. Developmental progression is governed by degree-days 
	 * (m_AgeDegrees) rather than chronological age for stages where temperature affects development.
	 */
	int m_Age;
	
	/**
	 * @var m_OurPopulationManager
	 * @brief Pointer to the Osmia population manager instance
	 * @details Provides access to population-level services: landscape queries, random number generation,
	 * object creation/deletion, and configuration parameters. Set at object creation and maintained
	 * through metamorphosis.
	 */
	Osmia_Population_Manager* m_OurPopulationManager;
	
	/**
	 * @var m_OurParasitoidPopulationManager
	 * @brief Static pointer to parasitoid population manager (if using mechanistic parasitoid model)
	 * @details Enables individuals to query parasitoid density when calculating parasitism risk for
	 * open nest cells. Static because there is only one parasitoid manager per simulation. May be
	 * nullptr if parasitism is calculated via simple probability model.
	 */
	static OsmiaParasitoid_Population_Manager* m_OurParasitoidPopulationManager;
	
	/**
	 * @var m_TempToday
	 * @brief Mean daily temperature (°C) for current timestep
	 * @details Updated once per day by population manager. Static because temperature is the same
	 * for all individuals at a given timestep. Used in degree-day calculations and temperature
	 * threshold tests.
	 * 
	 * @par Implementation Note
	 * Daily mean temperature is read from weather input file. Future versions may implement hourly
	 * temperatures for improved development accuracy under fluctuating conditions.
	 */
	static double m_TempToday;
	
	/**
	 * @var m_TempTodayInt
	 * @brief Rounded integer temperature for array indexing
	 * @details Computed as floor(m_TempToday + 0.5) for use in temperature-indexed lookup tables.
	 * Currently not extensively used, but available if temperature-dependent parameters are
	 * implemented as arrays rather than calculated values.
	 */
	static int m_TempTodayInt;
	
	/**
	 * @var m_DailyDevelopmentMortEggs
	 * @brief Daily mortality probability for eggs (probability per day)
	 * @details Default: 0.0014
	 * 
	 * @par Empirical Basis
	 * Based on Radmacher and Strohm (2011) observing 5.2% egg-to-cocoon mortality under fluctuating
	 * temperature regime (10-25°C), typical of field conditions. Value divided across developmental
	 * stages with equal allocation to egg and larva stages.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses the formal model value precisely. Formal model chose
	 * constant mortality despite observed temperature relationships because data were insufficient
	 * to reliably parameterize temperature-dependent functions.
	 * 
	 * @par Biological Interpretation
	 * Represents combined effects of: desiccation in hot/dry conditions, chilling injury at low
	 * temperatures, fungal/bacterial infection, and handling disturbance during laboratory studies.
	 * Field mortality may differ but is difficult to measure non-destructively.
	 * 
	 * @par Uncertainty
	 * HIGH - Laboratory studies may not capture full range of field mortality sources. Egg stage
	 * mortality shows high variation between studies (6-25%), possibly due to handling effects or
	 * environmental differences.
	 */
	static double m_DailyDevelopmentMortEggs;
	
	/**
	 * @var m_DailyDevelopmentMortLarvae
	 * @brief Daily mortality probability for larvae (probability per day)
	 * @details Default: 0.0014
	 * 
	 * @par Empirical Basis
	 * Same as egg mortality - derived from Radmacher and Strohm (2011) egg-to-cocoon value,
	 * with equal allocation to egg and larva stages. Giejdasz and Fliszkiewicz (2016) observed
	 * slightly higher larval mortality (12.6% total) but sample was smaller.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value without calibration.
	 * 
	 * @par Biological Interpretation
	 * Larval mortality primarily from: insufficient provision quality or quantity (Sedivy et al. 2011
	 * showed diet effects), fungal infection of provision mass, and potentially parasitoid larvae
	 * (tracked separately via m_ParasitoidStatus). Feeding larvae are relatively robust once
	 * established.
	 * 
	 * @par Uncertainty
	 * MEDIUM - More data available than for eggs, but field validation lacking. Provision quality
	 * effects not explicitly modelled beyond mass, potentially underestimating diet-related mortality.
	 */
	static double m_DailyDevelopmentMortLarvae;
	
	/**
	 * @var m_DailyDevelopmentMortPrepupae
	 * @brief Daily mortality probability for prepupae (probability per day)
	 * @details Default: 0.003
	 * 
	 * @par Empirical Basis
	 * Mean of Radmacher and Strohm (2011) and Giejdasz and Fliszkiewicz (2016) laboratory studies,
	 * both finding very low prepupal mortality (≤1.5%) across all temperature treatments. This
	 * slightly higher value provides conservative buffer for field uncertainty.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value.
	 * 
	 * @par Biological Interpretation
	 * Prepupae are diapausing, cocooned, and not feeding, making them relatively invulnerable to
	 * environmental stressors during the brief prepupal stage. Mortality primarily from pre-existing
	 * weakness (insufficient larval feeding) or cocoon failure allowing desiccation.
	 * 
	 * @par Uncertainty
	 * LOW - Consistent findings across multiple studies. This is the most reliably measured
	 * mortality parameter.
	 */
	static double m_DailyDevelopmentMortPrepupae;
	
	/**
	 * @var m_DailyDevelopmentMortPupae
	 * @brief Daily mortality probability for pupae (probability per day)
	 * @details Default: 0.003
	 * 
	 * @par Empirical Basis
	 * Identical to prepupal mortality - both laboratory studies found similarly low mortality
	 * (<1.5%) for both stages, with no clear difference between them.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value.
	 * 
	 * @par Biological Interpretation
	 * Like prepupae, pupae are cocooned and protected. Metamorphosis is energetically demanding but
	 * failures are rare if the larva was well provisioned. Most pupal mortality likely reflects
	 * developmental abnormalities rather than environmental stress.
	 * 
	 * @par Uncertainty
	 * LOW - Well-supported by laboratory data.
	 */
	static double m_DailyDevelopmentMortPupae;
	
	/**
	 * @var m_OsmiaEggDevelTotalDD
	 * @brief Total degree-days required for egg development to hatching
	 * @details Default: 86.0 degree-days above threshold
	 * 
	 * @par Empirical Basis
	 * Laboratory data suggested 37.0 DD with 13.8°C threshold (Giejdasz & Wilkaniec 2002), but
	 * implementation uses 86 DD with 0°C threshold for improved field realism.
	 * 
	 * @par Difference from Formal Model
	 * **MAJOR CALIBRATION** - Formal model specified 37.0 DD with LDT=13.8°C based on laboratory
	 * curve fitting. Implementation increased to 86 DD whilst lowering threshold to 0°C to achieve
	 * realistic field emergence timing. This compensatory adjustment maintains similar absolute
	 * development duration under typical spring temperatures whilst preventing unrealistic cessation
	 * of development at cooler field temperatures.
	 * 
	 * @par Biological Interpretation
	 * The high laboratory threshold (13.8°C) likely reflects experimental artifacts or overfitting
	 * to limited temperature ranges. Field-active *O. bicornis* readily develop at temperatures
	 * well below 13.8°C during early spring. The 0°C threshold is more biologically realistic,
	 * with compensatory increase in total DD maintaining appropriate timing.
	 * 
	 * @par Uncertainty
	 * MEDIUM - Calibration-derived value rather than direct measurement. Field validation of
	 * emergence timing supports the adjustment but direct degree-day observations in nests are lacking.
	 */
	static double m_OsmiaEggDevelTotalDD;
	
	/**
	 * @var m_OsmiaEggDevelThreshold
	 * @brief Lower developmental threshold temperature for eggs (°C)
	 * @details Default: 0.0°C
	 * 
	 * @par Empirical Basis
	 * Laboratory curve fitting suggested 13.8°C (Giejdasz & Wilkaniec 2002), but field observations
	 * indicate development proceeds at much cooler temperatures.
	 * 
	 * @par Difference from Formal Model
	 * **MAJOR CALIBRATION** - Formal model used 13.8°C from laboratory analysis. Implementation
	 * reduced to 0°C to allow development across full range of spring field temperatures. This
	 * prevents unrealistic developmental arrest during cool spring periods when bees are actively
	 * nesting.
	 * 
	 * @par Biological Rationale
	 * Zero represents a conservative biological minimum - true developmental cessation probably
	 * occurs slightly above freezing, but 0°C provides simple, robust threshold without requiring
	 * sub-zero temperature handling. Field data show *O. bicornis* successfully developing in nests
	 * experiencing temperatures as low as 5°C.
	 * 
	 * @par Uncertainty
	 * MEDIUM - Threshold choice interacts with total DD requirement, so validation must consider
	 * both parameters together. Emergence phenology data support the combined parameterisation.
	 */
	static double m_OsmiaEggDevelThreshold;
	
	/**
	 * @var m_OsmiaLarvaDevelTotalDD
	 * @brief Total degree-days required for larval development to prepupal stage
	 * @details Default: 422 degree-days above threshold
	 * 
	 * @par Empirical Basis
	 * Laboratory studies (Giejdasz & Wilkaniec 2002, Radmacher & Strohm 2011) suggested 422.4 DD
	 * with LDT=8.5°C. Implementation maintains the DD value but adjusts threshold.
	 * 
	 * @par Difference from Formal Model
	 * **MODERATE CALIBRATION** - Formal model: 422.4 DD with LDT=8.5°C. Implementation: 422 DD
	 * with LDT=4.5°C. The threshold reduction follows the same logic as egg parameters, allowing
	 * development at cooler field temperatures. Total DD is essentially unchanged (422 vs 422.4).
	 * 
	 * @par Biological Interpretation
	 * Larvae are feeding and growing rapidly, with high metabolic demands. Development rate responds
	 * strongly to temperature. The 4.5°C threshold is more consistent with field observations of
	 * larvae developing successfully during cool spring periods. Higher threshold would predict
	 * unrealistically long larval periods or failed development.
	 * 
	 * @par Uncertainty
	 * MEDIUM - Larval development shows less inter-study variation than egg stage, increasing
	 * confidence. However, provision quality effects (not explicitly modelled) may interact with
	 * temperature to affect actual development rates.
	 */
	static double m_OsmiaLarvaDevelTotalDD;
	
	/**
	 * @var m_OsmiaLarvaDevelThreshold
	 * @brief Lower developmental threshold temperature for larvae (°C)
	 * @details Default: 4.5°C
	 * 
	 * @par Difference from Formal Model
	 * **MODERATE CALIBRATION** - Reduced from 8.5°C to 4.5°C following same rationale as egg
	 * threshold adjustment. Laboratory-derived thresholds consistently overestimate field minima.
	 * 
	 * @par Biological Rationale
	 * The 4.5°C threshold better represents the temperature below which larval metabolic processes
	 * effectively cease. Feeding and digestion require active enzymatic processes that slow dramatically
	 * below this temperature but don't fully stop until near-freezing conditions.
	 */
	static double m_OsmiaLarvaDevelThreshold;
	
	/**
	 * @var m_OsmiaPupaDevelTotalDD
	 * @brief Total degree-days required for pupal development to adult eclosion
	 * @details Default: 570 degree-days above threshold
	 * 
	 * @par Empirical Basis
	 * Laboratory studies suggested 272.3 DD with LDT=13.2°C, but field calibration required major
	 * adjustment.
	 * 
	 * @par Difference from Formal Model
	 * **MAJOR CALIBRATION** - Formal model: 272.3 DD with LDT=13.2°C. Implementation: 570 DD with
	 * LDT=1.1°C. This represents the largest parameter adjustment in the model. Code comment notes:
	 * "changed from 13.2 to prevent pupal death" - indicating original parameters caused
	 * developmental failures under field temperature regimes.
	 * 
	 * @par Biological Rationale
	 * The dramatic increase in total DD compensates for much lower threshold, maintaining realistic
	 * absolute development duration. The 1.1°C threshold allows pupal development to proceed during
	 * cool summer periods that would otherwise cause developmental stalling with the 13.2°C threshold.
	 * Laboratory studies at constant temperatures may not capture the integration of development
	 * under naturally fluctuating conditions.
	 * 
	 * @par Implementation Note
	 * This calibration was essential for model functionality - original parameters led to widespread
	 * mortality because pupae couldn't accumulate sufficient DD under realistic summer temperature
	 * regimes in central Europe. The adjusted values produce emergence timing consistent with
	 * field observations.
	 * 
	 * @par Uncertainty
	 * MEDIUM-HIGH - This is a calibration-derived value with large departure from laboratory
	 * measurements. However, successful reproduction of field phenology validates the adjustment.
	 * More detailed nest temperature monitoring would improve parameterisation confidence.
	 */
	static double m_OsmiaPupaDevelTotalDD;
	
	/**
	 * @var m_OsmiaPupaDevelThreshold
	 * @brief Lower developmental threshold temperature for pupae (°C)
	 * @details Default: 1.1°C
	 * 
	 * @par Difference from Formal Model
	 * **MAJOR CALIBRATION** - Reduced from 13.2°C to 1.1°C. See m_OsmiaPupaDevelTotalDD for
	 * complete rationale - these two parameters were calibrated together.
	 */
	static double m_OsmiaPupaDevelThreshold;
	
	/**
	 * @var m_OsmiaPrepupalDevelTotalDays
	 * @brief Total days for prepupal development at optimal temperature
	 * @details Default: 45 days
	 * 
	 * @par Empirical Basis
	 * Laboratory studies show prepupal development is non-linear with temperature, with optimum
	 * around 22°C giving minimum ~24 days (Radmacher & Strohm 2011, Giejdasz & Fliszkiewicz 2016).
	 * 
	 * @par Difference from Formal Model
	 * **STRUCTURAL DIFFERENCE** - Formal model specified quadratic function with 24.3-day optimum
	 * at 22°C. Implementation uses simpler time-based approach with 45-day baseline and ±10%
	 * individual variation. This represents a fundamentally different developmental model structure.
	 * 
	 * @par Biological Rationale
	 * Prepupal diapause is complex, involving photoperiod independence and non-monotonic temperature
	 * response. The quadratic relationship is poorly constrained by available data and difficult to
	 * parameterize robustly. Time-based approach with temperature thresholds provides more stable
	 * model behaviour whilst capturing key biology: prepupae take ~1-2 months and respond to 
	 * temperature extremes but not in simple linear fashion.
	 * 
	 * @par Implementation Note
	 * The 45-day value represents nominal duration under moderate temperatures. Individual variation
	 * (±10%) creates spread in development times. Temperature affects development via threshold-based
	 * rules rather than rate modification: development proceeds above prewintering threshold (15°C)
	 * but is suspended below it.
	 * 
	 * @par Uncertainty
	 * HIGH - This is a pragmatic simplification of complex prepupal physiology. Future improvements
	 * could implement the formal model's quadratic function if additional data become available to
	 * robustly parameterise the non-linear response.
	 */
	static double m_OsmiaPrepupalDevelTotalDays;
	
	/**
	 * @var m_OsmiaPrepupalDevelTotalDays10pct
	 * @brief Pre-computed 10% of prepupal development time
	 * @details Computational efficiency variable storing m_OsmiaPrepupalDevelTotalDays * 0.1 to avoid
	 * repeated multiplication when applying individual variation (uniform distribution ±10% around
	 * nominal duration).
	 */
	static double m_OsmiaPrepupalDevelTotalDays10pct;
	
	/**
	 * @var m_OsmiaInCocoonOverwinteringTempThreshold
	 * @brief Temperature threshold (°C) for accumulating overwintering degree-days
	 * @details Default: 0.0°C
	 * 
	 * @par Biological Rationale
	 * During winter diapause proper, cocooned adults accumulate chilling at temperatures above
	 * freezing. This threshold defines when temperatures contribute to the chilling requirement
	 * needed for diapause completion.
	 * 
	 * @par Difference from Formal Model
	 * **IMPLEMENTATION DETAIL** - Formal model discussed three-phase overwintering but didn't
	 * specify all threshold values explicitly. Implementation provides operational values.
	 */
	static double m_OsmiaInCocoonOverwinteringTempThreshold;
	
	/**
	 * @var m_OsmiaInCocoonEmergenceTempThreshold
	 * @brief Temperature threshold (°C) for post-diapause emergence counter
	 * @details Default: 5.0°C (adjusted from original 12°C)
	 * 
	 * @par Biological Rationale
	 * After diapause completion, adults remain in cocoons until spring warming. This threshold
	 * determines when days count toward emergence trigger. Below this temperature, adults remain
	 * quiescent even if diapause is complete.
	 * 
	 * @par Implementation Note
	 * Code comment indicates original value was 12°C but was reduced to 5°C during calibration.
	 * Lower threshold allows earlier emergence in response to spring warming, better matching
	 * observed field phenology.
	 */
	static double m_OsmiaInCocoonEmergenceTempThreshold;
	
	/**
	 * @var m_OsmiaInCocoonPrewinteringTempThreshold
	 * @brief Temperature threshold (°C) for accumulating prewintering degree-days
	 * @details Default: 15.0°C
	 * 
	 * @par Empirical Basis
	 * Based on Sgolastra et al. (2011) using 15°C baseline for calculating prewinter DD accumulation
	 * in *O. lignaria*. Same baseline applied to *O. bicornis* in absence of species-specific data.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value precisely. This is the baseline
	 * temperature for equation: overwintering_mortality = 0.05 × DD_prewinter - 4.63
	 * 
	 * @par Biological Interpretation
	 * Temperatures above 15°C during late summer/autumn (prewintering period) keep prepupal metabolism
	 * elevated, depleting lipid reserves and reducing overwintering success. The 15°C threshold
	 * distinguishes warm (stressful) from cool (appropriate) prewinter conditions.
	 */
	static double m_OsmiaInCocoonPrewinteringTempThreshold;
	
	/**
	 * @var m_OsmiaInCocoonWinterMortConst
	 * @brief Intercept for overwintering mortality equation
	 * @details Default: -4.63
	 * 
	 * @par Empirical Basis
	 * From Sgolastra et al. (2011) linear regression relating *O. lignaria* male overwintering
	 * mortality to prewinter degree-day accumulation. Applied to both sexes in *O. bicornis*.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value without adjustment.
	 * 
	 * @par Biological Interpretation
	 * Equation: mortality_prob = 0.05 × DD_prewinter - 4.63. Negative intercept means low mortality
	 * at zero prewinter DD (cool, appropriate prewintering), with mortality increasing linearly as
	 * warm prewinter conditions accumulate degree-days.
	 */
	static double m_OsmiaInCocoonWinterMortConst;
	
	/**
	 * @var m_OsmiaInCocoonWinterMortSlope
	 * @brief Slope for overwintering mortality equation
	 * @details Default: 0.05
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Direct implementation of formal model specification.
	 * 
	 * @par Biological Interpretation
	 * Each degree-day of warm prewinter conditions increases mortality probability by 0.05. For
	 * example, 100 DD of warm prewintering gives: 0.05 × 100 - 4.63 = 0.37 mortality probability.
	 */
	static double m_OsmiaInCocoonWinterMortSlope;
	
	/**
	 * @var m_OsmiaInCocoonEmergCountConst
	 * @brief Intercept for emergence counter equation
	 * @details Default: 35.4819
	 * 
	 * @par Biological Function
	 * Part of equation determining when spring emergence occurs based on accumulated degree-days:
	 * emergence_counter = 35.4819 - 0.0147 × DD_accumulated. When counter reaches zero, adult
	 * emerges from nest.
	 * 
	 * @par Difference from Formal Model
	 * **NEAR MATCH** - Code comment shows original value 39.4819, adjusted to 35.4819. This minor
	 * calibration shifts emergence timing slightly earlier in spring.
	 */
	static double m_OsmiaInCocoonEmergCountConst;
	
	/**
	 * @var m_OsmiaInCocoonEmergCountSlope
	 * @brief Slope for emergence counter equation
	 * @details Default: -0.0147
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model value.
	 */
	static double m_OsmiaInCocoonEmergCountSlope;
	
	/**
	 * @var m_OsmiaFemaleMassFromProvMassConst
	 * @brief Intercept for calculating female mass from provision mass
	 * @details Default: 4.00 mg
	 * 
	 * @par Empirical Basis
	 * From Seidelmann (2010) empirical relationship for *O. bicornis*:
	 * female_mass = 0.25 × provision_mass + 4.00
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implementation uses formal model relationship precisely.
	 */
	static double m_OsmiaFemaleMassFromProvMassConst;
	
	/**
	 * @var m_OsmiaFemaleMassFromProvMassSlope
	 * @brief Slope for calculating female mass from provision mass
	 * @details Default: 0.25
	 * 
	 * @par Biological Interpretation
	 * Approximately 25% of provision mass is converted to bee biomass, with remaining 75% lost to
	 * metabolism, egestion, and cocoon construction. This conversion efficiency is consistent across
	 * provision mass range observed in *O. bicornis*.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Direct implementation of Seidelmann (2010) relationship.
	 */
	static double m_OsmiaFemaleMassFromProvMassSlope;

	/** @brief Minimum target provisioning mass for male nest cells (mg) */
	static double m_MaleMinTargetProvisionMass;
	
	/** @brief Maximum target provisioning mass for male nest cells (mg) */
	static double m_MaleMaxTargetProvisionMass;
	
	/** @brief Minimum target provisioning mass for female nest cells (mg) */
	static double m_FemaleMinTargetProvisionMass;
	
	/** @brief Maximum target provisioning mass for female nest cells (mg) */
	static double m_FemaleMaxTargetProvisionMass;
	
	/** @brief Maximum female adult mass (mg) - sets upper bound for provision mass calculations */
	static double m_FemaleMaxMass;
	
	/** @brief Minimum female adult mass (mg) - sets lower bound for provision mass calculations */
	static double m_FemaleMinMass;
	
	/** @brief Maximum male adult mass (mg) - males are smaller than females */
	static double m_MaleMaxMass;
	
	/**
	 * @var m_PollenScoreToMg
	 * @brief Conversion factor from pollen availability score to provisioned mass (mg pollen per day)
	 * @details Converts landscape-based pollen resource scores (arbitrary units reflecting flower
	 * density and quality) into actual provisioned pollen mass. Calibrated to produce realistic
	 * provisioning rates and nest completion times under typical landscape conditions.
	 */
	static double m_PollenScoreToMg;
	
	/**
	 * @var m_DensityDependentPollenRemovalConst
	 * @brief Parameter linking pollen depletion to *Osmia* population density
	 * @details Implements simple density dependence: as local *Osmia* density increases, pollen
	 * resources deplete faster through competition. Value determines strength of this effect.
	 * 
	 * @par Biological Rationale
	 * Multiple foraging *Osmia* females from the same area will compete for limited floral resources,
	 * reducing per-capita pollen collection rates at high densities. This provides negative feedback
	 * preventing unrealistic population growth.
	 */
	static double m_DensityDependentPollenRemovalConst;
	
	/** @brief Minimum time (days) required to construct and provision one nest cell */
	static double m_MinimumCellConstructionTime;
	
	/** @brief Maximum time (days) required to construct and provision one nest cell */
	static double m_MaximumCellConstructionTime;
	
	/**
	 * @var m_TotalNestsPossible
	 * @brief Maximum number of nests a female can complete in her lifetime
	 * @details Determines upper bound on reproductive output. Typical values 3-5 nests. Combined with
	 * eggs per nest (from Seidelmann 2010 relationship) determines lifetime fecundity potential.
	 * Actual nests completed depends on longevity, resource availability, and weather.
	 */
	static int m_TotalNestsPossible;
	
	/**
	 * @var m_BombylidProbability
	 * @brief Baseline probability of bombyliid fly parasitism per open nest cell
	 * @details Simple parasitism model variant: fixed probability per cell based on cell open duration.
	 * Alternative to mechanistic parasitoid population model.
	 */
	static double m_BombylidProbability;
	
	/**
	 * @var m_ParasitismProbToTimeCellOpen
	 * @brief Conversion factor relating cell open time (days) to parasitism probability
	 * @details Longer open cells have higher parasitism risk as they provide longer window for
	 * parasitoid discovery and attack. This parameter scales time to probability.
	 */
	static double m_ParasitismProbToTimeCellOpen;
	
	/**
	 * @var m_ParasitoidAttackChance
	 * @brief Per-capita parasitoid attack probabilities for mechanistic parasitoid model
	 * @details Vector holding attack chance parameters when using explicit parasitoid population
	 * dynamics. Size and structure depend on parasitoid population model implementation.
	 */
	static vector<double> m_ParasitoidAttackChance;
	
	/**
	 * @var m_OsmiaFemaleR50distance
	 * @brief Typical homing distance - distance at which 50% of females cannot return to nest (m)
	 * @details Default: 660m. Based on central place foraging literature for small bees. Used in
	 * movement probability distributions.
	 * 
	 * @par Biological Basis
	 * Derived from relationships between body size (intertegular span) and foraging range. Smaller
	 * bees have shorter effective foraging radii due to energetic constraints and navigation limits.
	 */
	static double m_OsmiaFemaleR50distance;
	
	/**
	 * @var m_OsmiaFemaleR90distance
	 * @brief Maximum homing distance - distance at which 90% of females cannot return (m)
	 * @details Default: 1430m. Represents extreme foraging range, used in dispersal movements and
	 * maximum resource search distances.
	 */
	static double m_OsmiaFemaleR90distance;
	
	/** @brief Duration of prenesting period after emergence (days) */
	static int m_OsmiaFemalePrenesting;
	
	/** @brief Maximum adult female lifespan (days) */
	static int m_OsmiaFemaleLifespan;
	
	/**
	 * @var m_generalmovementdistances
	 * @brief Probability distribution for foraging and nest-searching movements
	 * @details Pre-computed movement distance distribution (typically beta) matching R50/R90 parameters.
	 * Used for selecting movement distances during resource searches.
	 */
	static probability_distribution m_generalmovementdistances;
	
	/**
	 * @var m_dispersalmovementdistances
	 * @brief Probability distribution for dispersal movements (longer than foraging movements)
	 * @details Separate distribution for dispersal events when females seek new nesting areas. May
	 * have different shape than general movements to represent directed long-distance movements.
	 */
	static probability_distribution m_dispersalmovementdistances;
	
	/**
	 * @var m_eggspernestdistribution
	 * @brief Probability distribution for planned eggs per nest
	 * @details Generates stochastic variation in reproductive planning. Females "plan" egg number
	 * before beginning nest, then actual eggs laid may differ based on resource availability and
	 * mortality.
	 */
	static probability_distribution m_eggspernestdistribution;
	
	/**
	 * @var m_exp_ZeroToOne
	 * @brief Exponential-like probability distribution over range [0,1]
	 * @details Utility distribution for various stochastic processes requiring exponential-shaped
	 * probabilities over unit interval.
	 */
	static probability_distribution m_exp_ZeroToOne;
	
	/** @brief Mass conversion ratio from cocoon mass to provision mass required */
	static double m_CocoonToProvisionMass;
	
	/** @brief Mass conversion ratio from provision mass to resulting cocoon mass */
	static double m_ProvisionToCocoonMass;
	
	/** @brief Total provision mass loss from first to last egg in a nest (mg) */
	static double m_TotalProvisioningMassLoss;
	
	/** @brief Stochastic range around total provisioning mass loss (mg) */
	static double m_TotalProvisioningMassLossRange;
	
	/** @brief Pre-computed double of mass loss range for efficiency */
	static double m_TotalProvisioningMassLossRangeX2;
	
	/**
	 * @var m_UsingMechanisticParasitoids
	 * @brief Flag selecting parasitism model: true=mechanistic population model, false=simple probabilities
	 * @details Determines which parasitism calculation method is used. Mechanistic model tracks
	 * parasitoid populations explicitly; simple model uses fixed probabilities.
	 */
	static bool m_UsingMechanisticParasitoids;
	
	/**
	 * @var m_OsmiaFemaleBckMort
	 * @brief Daily background mortality for adult females outside nest (probability per day)
	 * @details Based on Giejdasz et al. (2016) finding 0.02 daily mortality under semi-natural
	 * conditions. Represents combined hazards of foraging, predation, weather exposure.
	 */
	static double m_OsmiaFemaleBckMort;
	
	/** @brief Minimum eggs planned per nest (sets lower bound for egg planning distribution) */
	static int m_OsmiaFemaleMinEggsPerNest;
	
	/** @brief Number of attempts allowed for finding suitable nest location before giving up */
	static int m_OsmiaFindNestAttemptNo;
	
	/** @brief Maximum eggs planned per nest (sets upper bound for egg planning distribution) */
	static int m_OsmiaFemaleMaxEggsPerNest;
	
	/**
	 * @var m_emergenceday
	 * @brief Probability distribution for day of emergence relative to population mean
	 * @details Creates stochastic spread in emergence dates across population. Based on field
	 * observations showing 2-3 week emergence period for *O. bicornis* populations.
	 */
	static probability_distribution m_emergenceday;
	
	/**
	 * @var m_ParasitoidStatus
	 * @brief Records parasitism status of this individual
	 * @details Set during egg/larval stages when parasitism events occur. Parasitised individuals
	 * die at prescribed time based on parasitoid type. Only one parasitoid type per individual.
	 */
	TTypeOfOsmiaParasitoids m_ParasitoidStatus;
	
	/**
	 * @var m_OurNest
	 * @brief Pointer to nest containing this individual
	 * @details Dual use: For stages egg through InCocoon, points to natal nest where individual
	 * develops. For adult females, points to nest currently being provisioned. Null when female
	 * is dispersing or searching for nest location.
	 */
	Osmia_Nest* m_OurNest;
	
	/**
	 * @var m_Mass
	 * @brief Mass of individual or provision (mg) - interpretation varies by life stage
	 * @details For eggs through pupae: provision mass in cell (determines adult size). For InCocoon
	 * and adults: adult body mass. For females during provisioning: current provision mass in cell
	 * under construction.
	 */
	double m_Mass;
	
	/**
	 * @var m_foragehours
	 * @brief Hours available for foraging today (decremented as foraging proceeds)
	 * @details Adult females have limited daily forage time due to weather, daylight, and other
	 * activities (nest construction, mating). This counter tracks remaining forage hours, preventing
	 * unrealistic within-day provisioning rates.
	 */
	int m_foragehours;

public:
	/**
	 * @brief Constructor creating new Osmia_Base individual
	 * @param data Pointer to struct_Osmia containing initialization data
	 * @details Initializes base attributes, sets population manager pointer, assigns initial state.
	 * Most specific initialization occurs in derived class constructors.
	 */
	Osmia_Base(struct_Osmia* data);
	
	/**
	 * @brief Reinitialize object from object pool
	 * @param data Pointer to struct_Osmia containing new initialization data
	 * @details Resets all attributes to initial state for object reuse. Object pooling avoids
	 * repeated allocation/deallocation overhead in large simulations with high turnover.
	 */
	void ReInit(struct_Osmia* data);
	
	/** @brief Virtual destructor ensuring proper cleanup in derived classes */
	virtual ~Osmia_Base();
	
	/**
	 * @brief Behavioural state for death - cleanup and removal
	 * @details Virtual method called when individual dies from mortality, parasitism, or old age.
	 * Handles any necessary cleanup before object deletion by population manager.
	 */
	virtual void st_Dying(void);
	
	/**
	 * @brief First phase of daily timestep
	 * @details Called once per day before Step loop. Empty in base class - derived classes override
	 * to implement stage-specific daily initialization (e.g., setting daily forage hours).
	 */
	virtual void BeginStep(void) { ; }
	
	/**
	 * @brief Main behavioural step - called repeatedly until DONE
	 * @details Core of state machine implementation. Empty in base class - derived classes implement
	 * specific behaviours. Population manager calls Step repeatedly until all individuals return
	 * terminal state (DONE or DIE).
	 */
	virtual void Step(void) { ; }
	
	/**
	 * @brief Final phase of daily timestep
	 * @details Called once per day after all Step calls complete. Empty in base class - derived
	 * classes use for daily wrap-up tasks (e.g., age incrementing, resource depletion).
	 */
	virtual void EndStep(void) { ; }
	
	/** @brief Get individual's age in days */
	int GetAge() { return m_Age; }
	
	/** @brief Set individual's age in days */
	void SetAge(int a_age) { m_Age = a_age; }
	
	/** @brief Get individual's mass (mg) - interpretation depends on life stage */
	double GetMass() { return m_Mass; }
	
	/** @brief Set individual's mass (mg) */
	void SetMass(double a_mass) { m_Mass = a_mass; }
	
	/**
	 * @brief Set parasitism status
	 * @param a_status Type of parasitoid affecting this individual
	 * @details Once set, parasitism status persists through metamorphosis. Parasitised individuals
	 * die at characteristic time for that parasitoid type.
	 */
	void SetParasitised(TTypeOfOsmiaParasitoids a_status) { 
		m_ParasitoidStatus = a_status; 
	}
	
	/** @brief Get parasitism status */
	TTypeOfOsmiaParasitoids GetParasitised(void) { return m_ParasitoidStatus; }
	
	/** @brief Get pointer to nest containing or being provisioned by this individual */
	Osmia_Nest* GetNest() { return m_OurNest; }
	
	/**
	 * @brief Populate all static parameters from configuration file
	 * @details Called once during population manager initialization. Reads configuration values
	 * and stores in static variables for efficient access during simulation. This avoids repeated
	 * configuration file lookups.
	 */
	static void SetParameterValues();
	
	/**
	 * @brief Set current daily temperature for all individuals
	 * @param a_temperature Mean daily temperature (°C)
	 * @details Static method called once per day by population manager. Updates m_TempToday and
	 * m_TempTodayInt for use in degree-day calculations across all individuals.
	 */
	static void SetTemp(double a_temperature) { 
		m_TempToday = a_temperature;
		m_TempTodayInt = int(floor(a_temperature + 0.5)); 
	}
	
	/**
	 * @brief Set parasitoid population manager pointer
	 * @param a_popman Pointer to parasitoid manager (or nullptr if not using mechanistic parasitoids)
	 * @details Static method called during initialization if mechanistic parasitoid model is active.
	 */
	static void SetParasitoidManager(OsmiaParasitoid_Population_Manager* a_popman) 
	{ 
		m_OurParasitoidPopulationManager = a_popman; 
	}
};

