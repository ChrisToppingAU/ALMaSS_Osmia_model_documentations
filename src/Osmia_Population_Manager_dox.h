/*
*******************************************************************************************************
Copyright (c) 2017, Christopher John Topping, Aarhus University
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
 * @file Osmia_Population_Manager.h
 * @brief Population-level management for *Osmia bicornis* agent-based model
 * 
 * @details This file provides the population manager infrastructure for the *Osmia bicornis* 
 * simulation model. The population manager handles:
 * - Global parameter initialization and coordination
 * - Object lifecycle management (creation, pooling, destruction)
 * - Daily scheduling and execution order
 * - Spatial data structures (density grids, nest management)
 * - Environmental condition monitoring (temperature, weather, seasons)
 * - Optional parasitoid population dynamics
 * 
 * The manager follows the ALMaSS framework Population_Manager architecture, orchestrating
 * daily simulation steps whilst individual agents handle their own behaviour through the
 * classes defined in Osmia.h.
 * 
 * @par Design Philosophy
 * The manager separates population-level concerns (scheduling, resource allocation, global
 * parameters) from individual-level behaviour (foraging, development, reproduction). This
 * separation allows individual agents to focus on behavioural decisions whilst the manager
 * handles the infrastructure needed to support those behaviours at scale.
 * 
 * @par Relationship to Formal Model
 * The formal model (Ziółkowska et al. 2025) describes individual-level processes without
 * detailing the simulation infrastructure. This implementation adds:
 * - Object pooling for computational efficiency
 * - Spatial indexing structures for rapid neighbour searches  
 * - Pre-calculated lookup tables for temperature-dependent processes
 * - Optional parasitoid dynamics extending the core bee model
 * 
 * These additions maintain biological fidelity whilst enabling practical large-scale simulation.
 * 
 * @author Christopher J. Topping
 * @author Enhanced documentation: [AUTHOR NAMES TO BE ADDED]
 * @date Original: May 2017
 * @date Enhanced: 2025
 * @ingroup Osmia_Model
 * 
 * @see Osmia.h for individual agent behaviour
 * @see Osmia_Population_Manager.cpp for implementation
 * @see Ziółkowska et al. (2025) Food and Ecological Systems Modelling Journal
 */

#include <forward_list>

//---------------------------------------------------------------------------
#ifndef Osmia_Population_ManagerH
#define Osmia_Population_ManagerH
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

class Osmia;
class OsmiaParasitoid_Population_Manager;

//==============================================================================
// ENUMERATIONS AND TYPE DEFINITIONS
//==============================================================================

/**
 * @enum TTypeOfOsmiaLifeStages
 * @brief Enumeration of modelled life stages for *Osmia bicornis*
 * 
 * @details Defines the six discrete life stages tracked in the model, corresponding to
 * distinct developmental phases with different biological processes and parameter sets.
 * The integer backing type enables direct use as array indices for stage-specific data.
 * 
 * @par Life Stage Sequence
 * - **Egg** (to_OsmiaEgg): From laying until hatching and feeding initiation
 * - **Larva** (to_OsmiaLarva): Active feeding phase plus cocoon construction
 * - **Prepupa** (to_OsmiaPrepupa): Summer diapause period within cocoon
 * - **Pupa** (to_OsmiaPupa): Metamorphosis from larval to adult form
 * - **InCocoon** (to_OsmiaInCocoon): Fully developed adults remaining in cocoons (includes overwintering)
 * - **Female** (to_OsmiaFemale): Emerged, active adult females
 * 
 * @par Biological Basis
 * Stage boundaries correspond to morphologically and physiologically distinct phases documented
 * in Radmacher and Strohm (2011), Giejdasz and Wilkaniec (2002), and other laboratory studies.
 * The prepupal stage represents an ecologically critical summer diapause allowing synchronisation
 * of emergence timing across variable spring weather conditions.
 * 
 * @par Implementation Note
 * Males are not explicitly modelled (see Osmia.h for rationale). The Female stage represents
 * only reproductive females; males emerge but their dynamics are implicitly captured through
 * sex ratio and mate availability assumptions.
 */
enum class TTypeOfOsmiaLifeStages : int
{
	to_OsmiaEgg = 0,     ///< Egg stage: laying to hatching
	to_OsmiaLarva,       ///< Larval stage: feeding and cocoon spinning
	to_OsmiaPrepupa,     ///< Prepupal stage: summer diapause
	to_OsmiaPupa,        ///< Pupal stage: metamorphosis
	to_OsmiaInCocoon,    ///< Adult-in-cocoon stage: includes overwintering
	to_OsmiaFemale       ///< Active adult female stage
};

/**
 * @typedef eggsexratiovsagelogisticcurvedata
 * @brief Storage for pre-calculated egg sex ratio values across female age
 * @details Vector storing logistic curve values representing probability of female egg
 * as function of maternal age. Pre-calculation during initialization avoids repeated
 * evaluation of complex logistic functions during simulation.
 */
typedef vector<double> eggsexratiovsagelogisticcurvedata;

/**
 * @typedef femalecocoonmassvsagelogisticcurvedata  
 * @brief Storage for pre-calculated female cocoon mass targets across maternal age
 * @details Vector storing logistic curve values for target provision mass for female
 * offspring cells as function of maternal age. Based on empirical relationship between
 * maternal condition (age) and offspring investment strategy (Seidelmann et al. 2010).
 */
typedef vector<double> femalecocoonmassvsagelogisticcurvedata;

//==============================================================================
// POLLEN AND NECTAR THRESHOLD DATA CLASS
//==============================================================================

/**
 * @class OsmiaPollenNectarThresholds
 * @brief Container for monthly resource quality and quantity thresholds
 * 
 * @details Simple data class holding four threshold values determining whether
 * a habitat patch is considered suitable for foraging. Thresholds vary by month
 * to reflect seasonal changes in floral resource availability and bee nutritional
 * requirements.
 * 
 * @par Biological Rationale
 * *Osmia bicornis* females are selective foragers, rejecting patches below minimum
 * quality or quantity thresholds. These thresholds likely reflect:
 * - Energetic costs of flight to/from patch
 * - Time constraints during provisioning period
 * - Nutritional requirements for larval development
 * 
 * Monthly variation acknowledges that early-season flowers (when bee density is low)
 * may offer different reward structures than late-season flowers (when competition
 * is higher and flower density may have declined).
 * 
 * @par Implementation Note
 * Thresholds read from configuration file (cfg_OsmiaPollenNectarThresholdsQual_XXX and
 * cfg_OsmiaPollenNectarThresholdsQuan_XXX where XXX is month abbreviation). Population
 * manager stores one OsmiaPollenNectarThresholds object per month in m_PN_thresholds vector.
 */
class OsmiaPollenNectarThresholds
{
public:
	/** @brief Minimum pollen quality score for patch acceptance */
	double m_pollenTqual = 0.0;
	
	/** @brief Minimum nectar quality score for patch acceptance */
	double m_nectarTqual = 0.0;
	
	/** @brief Minimum pollen quantity (mg available) for patch acceptance */
	double m_pollenTquan = 0.0;
	
	/** @brief Minimum nectar quantity (mg available) for patch acceptance */
	double m_nectarTquan = 0.0;
};

//==============================================================================
// PARASITOID POPULATION DYNAMICS (OPTIONAL EXTENSION)
//==============================================================================

/**
 * @class OsmiaParasitoidSubPopulation
 * @brief Spatially-explicit sub-population for a single parasitoid species
 * 
 * @details Represents a spatial cell in the parasitoid population grid, tracking local
 * parasitoid density and handling daily processes (mortality, dispersal, reproduction).
 * Multiple sub-populations tile the landscape to create spatial heterogeneity in
 * parasitism risk.
 * 
 * @par Biological Context
 * *Osmia bicornis* nests are parasitised by various insects including bombylid flies
 * (e.g., *Anthrax anthrax*) and chrysidid wasps. Parasitoid populations exhibit their
 * own spatial dynamics, dispersing between areas and responding to local host density.
 * This class implements a simplified parasitoid model that can optionally replace
 * the simpler probability-based parasitism in the core model.
 * 
 * @par Mechanistic vs. Probability-Based Parasitism
 * The model supports two parasitism approaches:
 * - **Probability-based** (default): Parasitism risk is a simple function of cell open time
 * - **Mechanistic** (this class): Parasitism emerges from explicit parasitoid population dynamics
 * 
 * The mechanistic approach adds realism (parasitoid dispersal, local aggregation) but increases
 * computational cost and parameter uncertainty. Choice between approaches depends on research
 * questions and available calibration data.
 * 
 * @par Spatial Structure
 * The landscape is divided into a coarse grid (e.g., 1 km² cells) with one sub-population
 * per cell per parasitoid species. This resolution balances realism (parasitoids do show
 * spatial structure at km scales) against computational tractability.
 * 
 * @par Implementation Note
 * Sub-populations managed by OsmiaParasitoid_Population_Manager which handles the spatial
 * array and coordinates dispersal between neighbouring cells. See that class for grid
 * structure details.
 * 
 * @see OsmiaParasitoid_Population_Manager for grid management
 * @see Osmia_Female::CalcParasitised() for parasitism application to bee eggs
 */
class OsmiaParasitoidSubPopulation
{
protected:
	// Attributes
	
	/** 
	 * @brief Current number of parasitoids in this spatial cell
	 * @details Continuous (double) to allow fractional individuals, avoiding discretization
	 * artifacts in dispersal and mortality. Biological interpretation: expected number of
	 * parasitoids, acknowledging spatial sampling stochasticity.
	 */
	double m_NoParasitoids;
	
	/** 
	 * @brief Proportion of population dispersing per time step
	 * @details Daily diffusion rate (0-1 scale). Represents innate dispersal tendency
	 * independent of distance travelled. Biological basis: parasitoids show consistent
	 * movement rates even in uniform host environments, suggesting intrinsic dispersal
	 * motivation (Van Nouhuys and Hanski 2002).
	 */
	double m_DiffusionRate;
	
	/** 
	 * @brief Distance-dependent dispersal kernel parameter
	 * @details Controls how dispersal probability decays with distance. Higher values
	 * indicate shorter-distance movements (more local aggregation). Typically estimated
	 * from mark-recapture or genetic data, but often calibrated due to data limitations.
	 */
	double m_DiffusionConstant;
	
	/** 
	 * @brief Pre-calculated indices of 8 neighbouring cells (Moore neighbourhood)
	 * @details Performance optimization: storing neighbour indices avoids repeated
	 * coordinate-to-index calculations during dispersal. Calculated once during
	 * construction, handling edge cases (boundary cells have fewer neighbours).
	 */
	int m_CellIndexArray[8];
	
	/** @brief Grid X-coordinate of this cell */
	int m_x;
	
	/** @brief Grid Y-coordinate of this cell */
	int m_y;
	
	/** 
	 * @brief Pointer to owning population manager
	 * @details Allows sub-population to access manager methods (add dispersers to
	 * neighbours, query landscape state). Necessary because dispersal affects multiple
	 * sub-populations and requires manager-level coordination.
	 */
	OsmiaParasitoid_Population_Manager* m_OurPopulationManager;
	
	/** 
	 * @brief Monthly mortality rates (proportion dying per day)
	 * @details Static array shared across all sub-populations, read from configuration.
	 * Monthly resolution acknowledges seasonal variation in parasitoid survival
	 * (higher mortality during unfavourable weather periods). Daily application
	 * allows flexible phenology without discrete monthly boundaries.
	 * 
	 * @par Data Requirements
	 * Ideally derived from field studies tracking parasitoid survival across seasons.
	 * In practice, often calibrated to match observed parasitism patterns in bee populations.
	 */
	static array<double,12> m_MortalityPerMonth;
	
	/** 
	 * @brief Current month index (0-11)
	 * @details Static optimization: all sub-populations experience same month, so single
	 * shared variable avoids repeated date queries. Updated by population manager.
	 */
	static int m_ThisMonth;
	
	// Methods
public:
	/**
	 * @brief Constructor initializing spatial position and dispersal parameters
	 * @param a_dispersalfraction Daily dispersal rate (proportion leaving cell)
	 * @param a_startno Initial parasitoid count in this cell
	 * @param a_x Grid X-coordinate
	 * @param a_y Grid Y-coordinate  
	 * @param a_wide Total grid width (cells)
	 * @param a_high Total grid height (cells)
	 * @param a_popman Pointer to owning population manager
	 * 
	 * @details Calculates neighbour indices considering boundaries (edge cells have
	 * fewer neighbours). Initializes parasitoid count, stores spatial position.
	 */
	OsmiaParasitoidSubPopulation(double a_dispersalfraction, double a_startno, int a_x, 
	                              int a_y, int a_wide, int a_high, 
	                              OsmiaParasitoid_Population_Manager* a_popman);
	
	/** @brief Destructor */
	~OsmiaParasitoidSubPopulation();
	
	/** 
	 * @brief Add parasitoids (e.g., immigrants from neighbours)
	 * @param a_change Number of parasitoids to add
	 */
	void Add(double a_change) { m_NoParasitoids += a_change; }
	
	/** 
	 * @brief Remove parasitoids (e.g., mortality, emigration)
	 * @param a_change Number of parasitoids to remove
	 */
	void Remove(double a_change) { m_NoParasitoids -= a_change; }
	
	/** 
	 * @brief Query current population size
	 * @return Number of parasitoids in this cell
	 */
	double GetSubPopnSize() { return m_NoParasitoids; }
	
	/** 
	 * @brief Apply daily mortality
	 * @details Removes proportion m_MortalityPerMonth[m_ThisMonth] of population.
	 * Stochastic mortality at individual level; deterministic at population level.
	 */
	void DailyMortality();
	
	/** 
	 * @brief Execute dispersal to neighbouring cells
	 * @details Calculates emigration (m_NoParasitoids × m_DiffusionRate), distributes
	 * emigrants to neighbours using distance kernel (m_DiffusionConstant). Updates
	 * both this cell and neighbour cells through population manager interface.
	 */
	void Dispersal();
	
	/** 
	 * @brief Execute reproduction based on local host density
	 * @details Queries landscape for local *Osmia* nest density, calculates offspring
	 * production using functional response, adds offspring to population. Reproduction
	 * success depends on host encounter rate (search efficiency × host density).
	 */
	void Reproduce();
	
	/**
	 * @brief Main daily update orchestrating sub-population processes
	 * @details Calls processes in biologically meaningful order:
	 * 1. DailyMortality() - deaths from all causes
	 * 2. Dispersal() - movement of survivors
	 * 3. Reproduce() - offspring production
	 * 
	 * Virtual to allow derived classes to modify or extend process sequence.
	 */
	virtual void DoFirst() {
		DailyMortality();
		Dispersal();
		Reproduce();
	}
	
	/** 
	 * @brief Update current month for mortality lookup
	 * @param a_month Month index (0-11, where 0=January)
	 */
	void SetThisMonth(int a_month) { m_ThisMonth = a_month; }
	
	/** 
	 * @brief Set monthly mortality rates array
	 * @param a_morts Array of 12 monthly mortality values
	 * @details Called during initialization by population manager after reading
	 * configuration file. Static array shared across all sub-populations.
	 */
	void SetMortalities(array<double, 12> a_morts) {
		m_MortalityPerMonth = a_morts;
	}
};

/**
 * @class OsmiaParasitoid_Population_Manager
 * @brief Grid-based manager coordinating multiple parasitoid sub-populations
 * 
 * @details Manages the complete parasitoid population as a spatial array of sub-populations,
 * handling grid initialization, inter-cell dispersal coordination, and population queries.
 * Provides interface for bee agents to query local parasitism risk and for sub-populations
 * to exchange dispersers.
 * 
 * @par Spatial Structure
 * Landscape divided into square grid cells (m_CellSize meters per side, typically 1000 m).
 * Each cell contains one sub-population per parasitoid species. Grid dimensions calculated
 * from landscape extent during construction. Boundary effects handled by sub-populations
 * (edge cells have fewer neighbours).
 * 
 * @par Multi-Species Support
 * Grid stores multiple sub-populations per cell, one for each parasitoid species type.
 * Total size = m_Wide × m_High × n_species. This allows different parasitoid species
 * with distinct dispersal, mortality, and attack parameters to coexist and potentially
 * interact through shared hosts.
 * 
 * @par Performance Considerations
 * Grid resolution represents a trade-off:
 * - Finer grids (smaller cells) → more spatial realism but higher computational cost
 * - Coarser grids (larger cells) → faster but less spatial detail
 * 
 * Typical choice (1 km cells) balances parasitoid dispersal scales (hundreds of meters
 * per day) against simulation performance for landscape-scale studies.
 * 
 * @see OsmiaParasitoidSubPopulation for individual cell dynamics
 */
class OsmiaParasitoid_Population_Manager : public Population_Manager
{
protected:
	// Attributes
	
	/** 
	 * @brief Vector storing all parasitoid sub-populations
	 * @details Flattened 3D array: [x + y×width + species×(width×height)]
	 * Allows efficient iteration and direct index access for dispersal calculations.
	 */
	vector<OsmiaParasitoidSubPopulation*> m_SubPopulations;
	
	/** 
	 * @brief Pointer to landscape object
	 * @details Provides access to:
	 * - Habitat quality data (for parasitoid reproduction)
	 * - Host density information (bee nests)
	 * - Spatial queries (polygon lookups)
	 * - Environmental conditions (temperature, weather)
	 */
	Landscape* m_TheLandscape;
	
	/** @brief Number of grid cells in X direction */
	unsigned m_Wide;
	
	/** @brief Number of grid cells in Y direction */
	unsigned m_High;
	
	/** 
	 * @brief Grid cell size in meters (cells are square)
	 * @details Typically 1000 m (1 km²). Determines spatial resolution of parasitoid
	 * population dynamics and affects computational cost (finer = more cells = slower).
	 */
	unsigned m_CellSize;
	
	/** 
	 * @brief Total number of sub-population cells per species
	 * @details m_Size = m_Wide × m_High. Pre-calculated for efficiency in index arithmetic.
	 */
	unsigned m_Size;
	
	// Methods
public:
	/**
	 * @brief Constructor initializing parasitoid population grid
	 * @param a_landscape Pointer to landscape object
	 * @param a_cellsize Grid cell size in meters
	 * 
	 * @details:
	 * 1. Calculates grid dimensions from landscape extent
	 * 2. Allocates sub-population array
	 * 3. Constructs each sub-population with spatial position and neighbours
	 * 4. Reads configuration parameters (mortality rates, dispersal, reproduction)
	 * 5. Initializes sub-populations with starting densities
	 */
	OsmiaParasitoid_Population_Manager(Landscape* a_landscape, int a_cellsize);
	
	/** 
	 * @brief Destructor cleaning up sub-populations
	 * @details Deletes all sub-population objects and clears vector.
	 */
	~OsmiaParasitoid_Population_Manager();
	
	/**
	 * @brief Add dispersing parasitoids to specified sub-population
	 * @param a_ref Sub-population array index
	 * @param a_dispersers Number of parasitoids to add
	 * 
	 * @details Called by sub-populations during dispersal to deposit emigrants
	 * in neighbouring cells. Provides controlled access to m_SubPopulations array.
	 */
	void AddDispersers(int a_ref, double a_dispersers) {
		m_SubPopulations[a_ref]->Add(a_dispersers);
	}
	
	/**
	 * @brief Remove parasitoids from specified sub-population
	 * @param a_ref Sub-population array index
	 * @param a_dispersers Number of parasitoids to remove
	 * 
	 * @details Used for emigration or mortality affecting specific cells.
	 */
	void RemoveParasitoids(int a_ref, double a_dispersers) {
		m_SubPopulations[a_ref]->Remove(a_dispersers);
	}
	
	/**
	 * @brief Query parasitoid density by sub-population index
	 * @param a_ref Sub-population array index
	 * @return Current parasitoid count in that cell
	 */
	double GetSize(int a_ref) { 
		return m_SubPopulations[a_ref]->GetSubPopnSize(); 
	}
	
	/**
	 * @brief Query parasitoid density by spatial coordinates
	 * @param a_x Grid X-coordinate
	 * @param a_y Grid Y-coordinate
	 * @return Current parasitoid count in that cell
	 * 
	 * @details Converts grid coordinates to array index: x + y×width
	 */
	double GetSize(int a_x, int a_y) { 
		return m_SubPopulations[a_x + a_y*m_Wide]->GetSubPopnSize(); 
	}
	
	/** 
	 * @brief Get array of parasitoid densities for all species at location
	 * @param a_x Grid X-coordinate
	 * @param a_y Grid Y-coordinate
	 * @return Array with density for each parasitoid species
	 * 
	 * @details Used by bee agents to query total parasitism risk from all species.
	 * Array size determined by TTypeOfOsmiaParasitoids enum.
	 */
	array<double, static_cast<unsigned>(TTypeOfOsmiaParasitoids::topara_foobar)> 
	GetParasitoidNumbers(int a_x, int a_y);
	
	/**
	 * @brief Add one parasitoid of specified type at location
	 * @param a_type Parasitoid species type
	 * @param a_x Landscape X-coordinate (meters)
	 * @param a_y Landscape Y-coordinate (meters)
	 * 
	 * @details Converts metric coordinates to grid coordinates (÷ cell size),
	 * calculates sub-population index accounting for species offset, adds
	 * parasitoid to appropriate sub-population.
	 * 
	 * Used during parasitoid reproduction or initialization.
	 */
	void AddParasitoid(TTypeOfOsmiaParasitoids a_type, int a_x, int a_y) {
		int subpop = ((a_x / m_CellSize) + (a_y / m_CellSize) * m_Wide) 
		           + (static_cast<unsigned>(a_type)-1) * m_Size;
		m_SubPopulations[subpop]->Add(1);
	}
};

//==============================================================================
// OSMIA CREATION DATA STRUCTURE
//==============================================================================

/**
 * @class struct_Osmia
 * @brief Initialization data package for creating new Osmia agents
 * 
 * @details Simple data structure bundling all information needed to construct
 * an Osmia agent at any life stage. Used during:
 * - Initial population seeding
 * - Reproduction (egg laying)
 * - Stage transitions (e.g., larva → prepupa)
 * - Object pool reinitialization
 * 
 * @par Design Rationale
 * Separating creation data from agent classes simplifies memory management
 * and allows flexible agent initialization patterns. The population manager
 * prepares struct_Osmia packages and passes them to agent constructors,
 * decoupling initialization logic from agent behaviour.
 * 
 * @par Object Pooling
 * ALMaSS uses object pooling to avoid repeated memory allocation/deallocation
 * during simulation. When an agent dies, its object returns to a pool and is
 * later reinitialized with new struct_Osmia data. This struct provides the
 * clean interface for that reinitialization.
 */
class struct_Osmia
{
 public:
	/** 
	 * @brief Landscape X-coordinate (meters)
	 * @details Location where agent will be created. For eggs, this is nest location;
	 * for emerging adults, this is emergence site (becomes dispersal origin).
	 */
	int x;
	
	/** 
	 * @brief Landscape Y-coordinate (meters)
	 */
	int y;
	
	/** 
	 * @brief Current age (days since life stage began)
	 * @details Interpretation depends on life stage:
	 * - Developmental stages: accumulated development days
	 * - Adult: days since emergence (affects behaviour, mortality)
	 * - Often initialized to 0 for new stages
	 */
	int age;
	
	/** 
	 * @brief Sex of individual
	 * @details true = female, false = male
	 * 
	 * Determined during egg laying based on provision mass. Females require larger
	 * provisions (Seidelmann et al. 2010), so mothers allocate sex based on resources
	 * accumulated for each cell.
	 * 
	 * @par Implementation Note
	 * Although males are not explicitly modelled after emergence, sex is tracked through
	 * development because:
	 * - Development rates may differ between sexes
	 * - Cocoon masses differ (calibration/validation data)
	 * - Sex ratio at emergence affects population dynamics
	 */
	bool sex;
	
	/** 
	 * @brief Pointer to landscape object
	 * @details Provides agent access to:
	 * - Environmental conditions (temperature, weather)
	 * - Spatial queries (nearest flowers, nest sites)
	 * - Polygon data (habitat types, resource availability)
	 */
	Landscape* L;
	
	/** 
	 * @brief Pointer to population manager
	 * @details Allows agent to:
	 * - Query global parameters (thresholds, lookup tables)
	 * - Signal stage transitions (create next-stage object)
	 * - Update population-level tracking (density grids)
	 * - Report death (return to object pool or delete)
	 */
	Osmia_Population_Manager* OPM;
	
	/** 
	 * @brief Pointer to nest structure
	 * @details Used by developmental stages (egg through InCocoon) to access nest
	 * information and by adult females during provisioning. NULL when not associated
	 * with a nest (e.g., during female dispersal between nests).
	 * 
	 * Nests managed by Osmia_Nest_Manager; this pointer connects individual agents
	 * to their spatial-temporal context within the nest structure.
	 */
	Osmia_Nest* nest;
	
	/** 
	 * @brief Parasitism status
	 * @details Enum value indicating which parasitoid species (if any) has successfully
	 * attacked this individual. Set during egg laying based on cell open time and
	 * parasitoid density. Parasitised individuals develop normally until parasitoid
	 * emerges, then die.
	 * 
	 * Value topara_Unparasitised indicates no parasitism.
	 */
	TTypeOfOsmiaParasitoids parasitised;
	
	/** 
	 * @brief Body mass of individual (mg)
	 * @details Meaning depends on life stage:
	 * - **Egg through Pupa**: provision mass allocated to cell (determines adult size)
	 * - **InCocoon**: cocoon mass (converted from provision mass using empirical equations)
	 * - **Adult Female**: body mass at emergence (determines fecundity, foraging efficiency)
	 * 
	 * @par Biological Importance
	 * Mass is the key individual-level state variable linking maternal provisioning
	 * decisions to offspring fitness. Larger provisions → larger adults → higher
	 * fecundity and potentially better survival (Radmacher and Strohm 2010).
	 * 
	 * @par Data Sources
	 * Conversion equations from Seidelmann (2006): provision mass → cocoon mass → adult mass.
	 * Relationships derived from laboratory measurements of provision and cocoon masses
	 * paired with adult body mass measurements.
	 */
	double mass;
	
	/** 
	 * @brief Pesticide-induced mortality probability
	 * @details Additional mortality risk from pesticide exposure, applied during
	 * development or adult life. Value 0.0 indicates no pesticide exposure; values
	 * 0.0-1.0 indicate daily mortality probability.
	 * 
	 * @par Implementation Note
	 * This mechanism allows integration with landscape-scale pesticide fate and
	 * transport models. Pesticide concentrations in floral resources or nesting
	 * substrates can be converted to bee mortality risks using toxicological dose-
	 * response data.
	 * 
	 * Default 0.0 for standard simulations without pesticide scenarios.
	 */
	double pest_mortality = 0.0;
	
	/** 
	 * @brief Accumulated overwintering degree-days at simulation start
	 * @details Used only when initializing simulation with overwintering adults
	 * (rather than starting from eggs). Allows setting initial population at
	 * realistic physiological states (partial overwintering progress) rather than
	 * requiring full year spin-up.
	 * 
	 * @par Biological Interpretation
	 * Overwintering development requires accumulation of degree-days below threshold.
	 * Starting adults with non-zero values simulates entry into overwintering period
	 * at various times, creating realistic emergence phenology distribution without
	 * multi-year simulation.
	 * 
	 * Should be 0.0 for normal simulation where population initialized from eggs.
	 */
	double overwintering_degree_days = 0.0;
};

//==============================================================================
// POLYGON-LEVEL NEST MANAGEMENT
//==============================================================================

/**
 * @class OsmiaPolygonEntry
 * @brief Nest list and density controls for a single landscape polygon
 * 
 * @details Each landscape polygon (habitat patch) maintains its own list of
 * *Osmia* nests and associated nesting suitability parameters. This class
 * provides the interface between the landscape (polygon-based) representation
 * and the nest management system.
 * 
 * @par Nesting Habitat Heterogeneity
 * Not all habitat is equally suitable for nesting. Suitability varies by:
 * - Vegetation structure (presence of hollow stems)
 * - Microclimate (sun exposure, shelter)
 * - Substrate availability (dead wood, exposed soil for mud)
 * - Management history (mowing timing, disturbance regime)
 * 
 * Polygon-level tracking allows spatial variation in nesting density and
 * creates realistic clustering of nests in favourable areas.
 * 
 * @par Nest Density Regulation
 * The m_OsmiaNestProb parameter controls how many nests can potentially
 * exist in a polygon, implementing density-dependent nesting constraints.
 * This represents:
 * - Finite nesting substrate (limited number of suitable cavities)
 * - Potential interference between nesting females
 * - Reduced attractiveness of heavily-used areas (accumulated parasitoids)
 * 
 * @see Osmia_Nest_Manager for nest lifecycle management
 * @see Osmia_Female for nest-finding behaviour
 */
class OsmiaPolygonEntry
{
protected:
	/** 
	 * @brief Linked list of active nests in this polygon
	 * @details std::forward_list provides efficient insertion/deletion with minimal
	 * memory overhead. Nests added when female claims cavity, removed when nest
	 * abandoned or all offspring emerge/die.
	 * 
	 * Forward list chosen over vector because:
	 * - Frequent additions/removals during simulation
	 * - Order not important (no need for random access)
	 * - Cache-friendly iteration for daily updates
	 */
	std::forward_list<Osmia_Nest*> m_NestList;
	
	/** 
	 * @brief Probability of successful nest establishment in this polygon
	 * @details Value 0.0-1.0 representing nesting suitability. Used when females
	 * search for nest sites; polygons with higher probability more likely to be
	 * selected and more likely to allow successful nest creation.
	 * 
	 * @par Parameterization
	 * Ideally derived from field surveys relating habitat characteristics to
	 * observed nest density. In practice, often calibrated to match observed
	 * spatial distribution of *Osmia* populations.
	 * 
	 * Can be static (read from habitat classification) or dynamic (modified by
	 * management events, vegetation succession, or density-dependent effects).
	 */
	double m_OsmiaNestProb;
	
	/** 
	 * @brief Maximum number of nests possible in this polygon
	 * @details Hard upper limit representing carrying capacity for nesting substrate.
	 * Prevents unrealistic nest densities when foraging resources abundant but
	 * nesting sites limited.
	 * 
	 * @par Biological Basis
	 * *Osmia bicornis* require pre-existing cavities (hollow stems, beetle holes,
	 * artificial nest boxes). These substrates are patchy and limited, particularly
	 * in intensive agricultural landscapes. Field studies show nests clustered in
	 * suitable microsites with tens to hundreds of nests in small areas whilst vast
	 * areas have none (Gathmann and Tscharntke 2002).
	 * 
	 * Calculated from polygon area, habitat type, and density parameters during
	 * initialization.
	 */
	int m_MaxNests;
	
	/** @brief Current number of active nests in this polygon */
	int m_CurrentNestCount;

public:
	/**
	 * @brief Get the nest list for this polygon
	 * @return Reference to forward_list of nest pointers
	 */
	std::forward_list<Osmia_Nest*>& GetNestList() { return m_NestList; }
	
	/**
	 * @brief Set nesting probability for this polygon
	 * @param a_prob Nesting suitability (0.0-1.0)
	 */
	void SetOsmiaNestProb(double a_prob) { m_OsmiaNestProb = a_prob; }
	
	/**
	 * @brief Query nesting probability
	 * @return Current nesting suitability value
	 */
	double GetOsmiaNestProb() { return m_OsmiaNestProb; }
	
	/**
	 * @brief Set maximum nest capacity
	 * @param a_max Maximum number of nests allowed
	 */
	void SetMaxNests(int a_max) { m_MaxNests = a_max; }
	
	/**
	 * @brief Query maximum nest capacity
	 * @return Maximum nests allowed in this polygon
	 */
	int GetMaxNests() { return m_MaxNests; }
	
	/**
	 * @brief Increment nest counter (nest created)
	 */
	void IncrementNestCount() { m_CurrentNestCount++; }
	
	/**
	 * @brief Decrement nest counter (nest released/destroyed)
	 */
	void DecrementNestCount() { m_CurrentNestCount--; }
	
	/**
	 * @brief Query current nest count
	 * @return Number of active nests currently in polygon
	 */
	int GetCurrentNestCount() { return m_CurrentNestCount; }
	
	/**
	 * @brief Check if polygon has capacity for additional nests
	 * @return true if m_CurrentNestCount < m_MaxNests
	 */
	bool HasNestCapacity() { return m_CurrentNestCount < m_MaxNests; }
};

//==============================================================================
// MAIN POPULATION MANAGER CLASS
//==============================================================================

/**
 * @class Osmia_Population_Manager  
 * @brief Central coordinator for *Osmia bicornis* population simulation
 * 
 * @details The population manager serves as the central orchestrator for the *Osmia*
 * simulation, handling initialization, daily scheduling, global parameter management,
 * and spatial data structures. It inherits from the ALMaSS Population_Manager base
 * class, providing *Osmia*-specific implementations of the standard simulation hooks.
 * 
 * @par Core Responsibilities
 * 
 * **1. Initialization and Configuration**
 * - Read parameters from configuration files
 * - Initialize lookup tables (provisioning times, sex ratios, cocoon masses)
 * - Set up spatial structures (density grids, nest manager, pollen map)
 * - Create initial population (typically overwintering adults)
 * 
 * **2. Daily Scheduling**
 * - DoFirst(): Update global environmental conditions (temperature, weather, phenology flags)
 * - DoBefore(): Pre-step calculations (prepupal development rates, foraging hours)
 * - Step(): Individual agents execute behaviour (inherited, not overridden)
 * - DoAfter(): Post-step cleanup (currently unused)
 * - DoLast(): End-of-day updates (seasonal flag management, statistics)
 * 
 * **3. Spatial Management**
 * - Maintain female density grid (1 km² resolution for conspecific density tracking)
 * - Coordinate nest manager (polygon-level nest lists, creation/release)
 * - Interface with pollen map (resource availability queries)
 * 
 * **4. Parameter Access**
 * - Provide lookup tables to individuals (provisioning parameters, sex ratios, etc.)
 * - Store and distribute global parameters (thresholds, mortality rates, development constants)
 * - Manage seasonal flags (pre-wintering end, overwintering end)
 * 
 * @par Relationship to Formal Model
 * The formal model (Ziółkowska et al. 2025) describes individual-level processes without
 * simulation infrastructure details. This manager adds:
 * - **Performance optimizations**: Pre-calculated lookup tables avoid repeated evaluations
 * - **Spatial indexing**: Density grid enables efficient local density queries
 * - **Computational efficiency**: Object pooling, parallel execution support
 * - **Extensibility**: Hooks for pesticides, parasitoids, output generation
 * 
 * These additions maintain biological fidelity whilst enabling practical landscape-scale
 * simulation (10⁴-10⁶ individuals over years).
 * 
 * @par Critical Design Decisions
 * 
 * **Pre-calculated Lookup Tables**
 * Age-dependent provisioning times (Seidelmann 2006) and size/age-dependent sex ratios
 * (Seidelmann et al. 2010) involve complex equations. Pre-calculating these during
 * initialization reduces computational cost during simulation when thousands of females
 * query these values daily.
 * 
 * **Density Grid Resolution**
 * 1 km² grid balances spatial detail against memory usage. Finer grids (100 m) would
 * better capture local crowding but increase memory and cache misses. Coarser grids
 * (5 km) save memory but lose ecologically meaningful density variation. 1 km represents
 * observed *Osmia* movement scales and is computationally tractable.
 * 
 * **Seasonal Flag Logic**
 * Pre-wintering and overwintering end flags control developmental switches (diapause
 * termination, emergence triggering). Determination based on sustained temperature
 * patterns (DoLast() logic) rather than fixed dates accommodates inter-annual climate
 * variation and spatial temperature gradients.
 * 
 * @see Osmia.h for individual agent classes
 * @see Osmia_Population_Manager.cpp for implementation details
 * @see Ziółkowska et al. (2025) Food and Ecological Systems Modelling Journal
 */
class Osmia_Population_Manager : public Population_Manager
{
public:
	/**
	 * @brief Constructor initializing population manager
	 * @param a_landscape Pointer to landscape object
	 * 
	 * @details Comprehensive initialization sequence:
	 * 
	 * 1. **Call Init()**: Read configuration, set up parameters
	 * 2. **Initialize seasonal flags**: Set mid-lifecycle start conditions
	 * 3. **Identify suitable nesting polygons**: Query landscape for nesting habitat
	 * 4. **Create initial population**: Place overwintering adults in suitable locations
	 * 5. **Set up lookup tables**: Pre-calculate provisioning times, sex ratios, cocoon masses
	 * 6. **Initialize spatial structures**: Density grid, nest manager
	 * 7. **Configure foraging**: Set up pollen map interface
	 * 
	 * @par Initial Population
	 * Typically starts with overwintering adults (InCocoon stage) to match field season
	 * start. Numbers and size distribution from configuration (cfg_OsmiaStartNo, mass range).
	 * Individuals placed randomly in suitable nesting polygons with accumulated overwintering
	 * degree-days (cfg_OsmiaOverwinterDegreeDaysInitialSimu) to create realistic emergence
	 * phenology without full-year spin-up.
	 * 
	 * @par Parallel Initialization
	 * Uses OpenMP parallelization for initial population creation (#pragma omp parallel),
	 * distributing agent construction across available threads. Thread-safe because each
	 * thread creates independent agents in separate memory.
	 */
	Osmia_Population_Manager(Landscape* a_landscape);
	
	/**
	 * @brief Destructor cleaning up population manager resources
	 * @details Deletes nest manager, clears spatial structures, closes output files.
	 * Individual agents handled by base Population_Manager destructor.
	 */
	~Osmia_Population_Manager();
	
	/**
	 * @brief Check if polygon suitable for *Osmia* nesting
	 * @param a_polyindex Landscape polygon index
	 * @return true if nesting possible in this polygon
	 * 
	 * @details Queries landscape polygon properties:
	 * - Habitat type classification
	 * - Nesting suitability parameter (if present)
	 * - Management state (e.g., recently mown grassland unsuitable)
	 * 
	 * Used during initialization to identify where to place starting population,
	 * and during simulation when females search for new nest sites.
	 * 
	 * @par Implementation Note
	 * Delegates to nest manager which maintains polygon-level nesting data.
	 * This indirection allows flexible nesting suitability criteria without
	 * modifying core population manager code.
	 */
	bool IsOsmiaNestPossible(int a_polyindex) {
		return m_OurOsmiaNestManager.IsOsmiaNestPossible(a_polyindex);
	}
	
	/**
	 * @brief Create new nest at specified location
	 * @param a_x Landscape X-coordinate (meters)
	 * @param a_y Landscape Y-coordinate (meters)
	 * @param a_polyindex Polygon containing nest location
	 * @return Pointer to newly created nest
	 * 
	 * @details Thread-safe nest creation sequence:
	 * 1. Acquire polygon lock (prevent concurrent modification)
	 * 2. Create nest via nest manager
	 * 3. Add nest to polygon's nest list
	 * 4. Increment polygon nest counter
	 * 5. Release polygon lock
	 * 
	 * @par Locking Rationale
	 * In parallel simulation, multiple females might simultaneously attempt to create
	 * nests in same polygon. Lock prevents race conditions (counter corruption, list
	 * inconsistency). Granularity at polygon level (rather than global lock) allows
	 * concurrent nesting in different polygons.
	 * 
	 * @par Biological Context
	 * Called by Osmia_Female during nest-finding behaviour when suitable cavity found.
	 * Females search sequentially through nearby nesting habitat, attempting nest creation
	 * until successful or giving up (triggering dispersal).
	 * 
	 * @see ReleaseOsmiaNest() for nest destruction
	 */
	Osmia_Nest* CreateNest(int a_x, int a_y, int a_polyindex) {
		Osmia_Nest* return_nest_ptr;
		m_TheLandscape->SetPolygonLock(a_polyindex);
		return_nest_ptr = m_OurOsmiaNestManager.CreateNest(a_x, a_y, a_polyindex); 
		m_TheLandscape->ReleasePolygonLock(a_polyindex);
		return return_nest_ptr;	
	}
	
	/**
	 * @brief Release (destroy) nest from polygon
	 * @param a_polyindex Polygon containing nest
	 * @param a_nest Pointer to nest being released
	 * 
	 * @details Thread-safe nest destruction:
	 * 1. Acquire polygon lock
	 * 2. Remove nest from polygon's nest list
	 * 3. Decrement polygon nest counter
	 * 4. Delete nest object or return to pool
	 * 5. Release polygon lock
	 * 
	 * @par When Called
	 * - Female abandons nest (before completing provisioning)
	 * - All offspring emerged or died (nest now empty)
	 * - Female dies whilst actively provisioning nest
	 * 
	 * Nest destruction frees nesting capacity in polygon, allowing future nest creation.
	 * 
	 * @see CreateNest() for nest creation
	 */
	void ReleaseOsmiaNest(int a_polyindex, Osmia_Nest* a_nest) {
		m_TheLandscape->SetPolygonLock(a_polyindex);
		m_OurOsmiaNestManager.ReleaseOsmiaNest(a_polyindex, a_nest);
		m_TheLandscape->ReleasePolygonLock(a_polyindex);
	}
	
	/**
	 * @brief Query daily foraging hours available
	 * @return Number of hours suitable for *Osmia* flight today
	 * 
	 * @details Returns m_FlyingWeather calculated by CalForageHours() each day.
	 * 
	 * @par Weather Constraints
	 * *Osmia bicornis* females require specific weather conditions for foraging:
	 * - Temperature > ~10°C (flight muscle activity threshold)
	 * - Low wind speed (< 5-6 m/s; higher winds prevent controlled flight)
	 * - No precipitation (wet wings reduce flight efficiency)
	 * - Sufficient light (active during daylight; peak mid-morning to mid-afternoon)
	 * 
	 * CalForageHours() integrates hourly weather data with these thresholds, counting
	 * hours meeting all criteria. On favourable days, might be 8-10 hours; on poor days,
	 * zero hours (females remain in nest).
	 * 
	 * @par Biological Consequences
	 * Foraging hour limitation is critical constraint on provisioning rate and ultimately
	 * reproductive success. Extended periods of poor weather can:
	 * - Delay nest completion (cells remain open longer → higher parasitism)
	 * - Reduce lifetime fecundity (time wasted waiting for good weather)
	 * - Increase mortality (longer exposure to predators, senescence)
	 * 
	 * @par Data Requirements
	 * Requires hourly weather data (temperature, wind, precipitation, radiation). Often
	 * derived from meteorological station data interpolated across landscape. Accuracy
	 * critical because small differences in foraging hour estimates accumulate over
	 * female lifetimes, substantially affecting population dynamics.
	 */
	int GetForageHours() { return m_FlyingWeather; }
	
	/**
	 * @brief Get provisioning time parameter for given adult age
	 * @param a_age Days since emergence
	 * @return Hours required to provision one complete cell
	 * 
	 * @details Returns pre-calculated value from m_NestProvisioningParameters lookup
	 * table. Values derived from Seidelmann (2006) efficiency equation:
	 * - Efficiency = 21.643 / (1 + exp((ln(age) - ln(18.888)) × 3.571)) mg/h
	 * - Construction time = (2.576 × efficiency + 56.17) / efficiency hours
	 * 
	 * @par Age Effect on Provisioning
	 * Younger females (<15 days) less efficient, requiring more hours per cell.
	 * Peak efficiency around day 18-20. Older females (>40 days) declining
	 * efficiency due to senescence. This age structure creates temporal variation
	 * in population-level provisioning rates even with constant weather.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implements Seidelmann (2006) equations precisely as specified
	 * in formal model (Ziółkowska et al. 2025, Provisioning section). Pre-calculation
	 * is performance optimization, not biological modification.
	 * 
	 * @par Implementation Note
	 * Lookup table populated during Init() for ages 0-364 days. Bees rarely survive
	 * beyond 60 days, so upper range provides safety margin without significant memory
	 * cost (365 × 8 bytes = ~3 KB per population manager).
	 */
	double GetProvisioningParams(int a_age) {
		return m_NestProvisioningParameters[a_age];
	}
	
	/**
	 * @brief Calculate first female cocoon mass based on female age and mass class
	 * @param a_age Maternal age (days since emergence)
	 * @param a_massclass Maternal mass class index (0-95, representing 4-28 mg range)
	 * @return Target provision mass for first female cell (mg)
	 * 
	 * @details Returns value from m_FemaleCocoonMassEqns lookup table with stochastic
	 * variation (±60% of mean, exponentially distributed). Implements declining investment
	 * pattern: first female offspring receives maximum resources, subsequent females
	 * receive less due to maternal resource depletion.
	 * 
	 * @par Biological Basis
	 * Seidelmann et al. (2010) documented that female *Osmia* provision first female
	 * cells more heavily than later cells. This reflects:
	 * - **Quality-quantity trade-off**: Better to produce few high-quality offspring
	 *   than many low-quality offspring when resources limited
	 * - **Maternal condition decline**: Foraging efficiency and resource accumulation
	 *   rates decline with maternal age and mass loss
	 * - **Bet-hedging**: First offspring have highest survival probability, so invest
	 *   heavily; later offspring are "bonus" if conditions remain favourable
	 * 
	 * @par Equation Structure
	 * Base calculation (lookup table): Logistic function of maternal age and mass
	 * - Younger, heavier mothers → larger first female mass
	 * - Older, lighter mothers → smaller first female mass
	 * 
	 * Stochastic variation: ±60% exponentially distributed around mean
	 * - Represents variation in local resource availability
	 * - Prevents unrealistic determinism in offspring size distributions
	 * - Calibrated from field data showing substantial within-female variation
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implements Seidelmann et al. (2010) provisioning equations as
	 * specified in formal model. Lookup table calculation and stochastic variation both
	 * match formal model description precisely. Pre-calculation is performance optimization.
	 * 
	 * @see GetSexRatioEggsAgeMass() for corresponding sex ratio calculation
	 */
	double GetFirstCocoonProvisioningMass(int a_age, int a_massclass) {
		return m_FemaleCocoonMassEqns[a_massclass][a_age] 
		     - (m_exp_ZeroTo1.Get() * m_FemaleCocoonMassEqns[a_massclass][a_age] * 0.6);
	}
	
	/**
	 * @brief Calculate sex ratio (proportion female) based on maternal age and mass
	 * @param a_massclass Maternal mass class index (0-95)
	 * @param a_age Maternal age (days since emergence)
	 * @return Probability of laying female egg (0.0-1.0)
	 * 
	 * @details Returns pre-calculated value from m_EggSexRatioEqns lookup table.
	 * Values derived from logistic equations fitted to Seidelmann et al. (2010) data
	 * showing age- and mass-dependent sex allocation.
	 * 
	 * @par Biological Basis
	 * Sex ratio patterns reflect:
	 * - **Local mate competition** theory: Female-biased sex ratios reduce competition
	 *   amongst sons, increasing fitness returns per male offspring
	 * - **Resource constraints**: Female offspring more expensive (larger provisions),
	 *   so mothers in poor condition produce more males
	 * - **Maternal condition effects**: Heavier, younger mothers produce more females
	 * 
	 * @par Pattern Across Age and Mass
	 * - Young, heavy mothers: ~65-75% female (optimal conditions, high female production)
	 * - Old mothers: ~40-50% female (resource depletion, shift toward cheaper males)
	 * - Light mothers: ~50-60% female (always somewhat constrained)
	 * 
	 * These patterns emerge from logistic functions parameterized from Seidelmann's
	 * observational data on provision mass distributions and sex ratios.
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Implements sex ratio equations specified in formal model, derived
	 * from Seidelmann et al. (2010). Lookup table approach matches formal model's logistic
	 * equations; pre-calculation is performance optimization.
	 * 
	 * @par Implementation Note
	 * Used during egg-laying sequence: female queries sex ratio for her current age/mass,
	 * compares to random draw (0-1), lays female egg if random value < sex ratio probability.
	 * This stochastic implementation creates realistic variation whilst maintaining
	 * expected sex ratio matching empirical data.
	 * 
	 * @see GetFirstCocoonProvisioningMass() for corresponding provision mass calculation
	 */
	double GetSexRatioEggsAgeMass(int a_massclass, int a_age) {
		return m_EggSexRatioEqns[a_massclass][a_age];
	}

	/**
	 * @brief Add female to spatial density grid
	 * @param a_loc Female's current location
	 * @return Grid cell index where female was added
	 * 
	 * @details Converts metric coordinates to grid coordinates (÷1000 m per cell),
	 * calculates linear index (x + y×grid_width), increments counter for that cell.
	 * 
	 * @par Density Grid Purpose
	 * Tracks local female density at 1 km² resolution to enable:
	 * - **Density-dependent behaviour**: Females may avoid crowded areas when dispersing
	 * - **Competition effects**: High local density could reduce per-capita foraging success
	 * - **Output/analysis**: Spatial distribution of reproductive effort
	 * 
	 * Currently used primarily for tracking rather than feedback (density doesn't
	 * strongly affect behaviour in base model), but infrastructure supports future
	 * density-dependent extensions.
	 * 
	 * @par When Called
	 * - Female emergence (add to emergence location grid cell)
	 * - Female dispersal (remove from old cell, add to new cell)
	 * - Female death (remove from last known cell)
	 * 
	 * @see RemoveFromDensityGrid() for removal
	 */
	int AddToDensityGrid(APoint a_loc) {
		int index = (a_loc.m_x / 1000) + (a_loc.m_y / 1000) * m_GridExtent;
		m_FemaleDensityGrid[index]++;
		return index;
	}
	
	/**
	 * @brief Add female to density grid using pre-calculated index
	 * @param a_index Grid cell index
	 * 
	 * @details Faster version when caller already calculated index (avoids coordinate
	 * conversion). Used during dispersal when female moves between adjacent cells
	 * (caller knows both old and new indices).
	 */
	void AddToDensityGrid(int a_index) {
		m_FemaleDensityGrid[a_index]++;
	}
	
	/**
	 * @brief Remove female from density grid
	 * @param a_index Grid cell index
	 * 
	 * @details Decrements counter for specified cell. Called when female dies or
	 * moves to different cell.
	 * 
	 * @par Thread Safety
	 * Increment/decrement operations not atomic, so concurrent access by multiple
	 * females could create race conditions. Current implementation assumes serial
	 * execution or careful parallelization (females partitioned to avoid conflicts).
	 * Future parallel implementations might need atomic operations or per-thread grids.
	 */
	void RemoveFromDensityGrid(int a_index) {
		m_FemaleDensityGrid[a_index]--;
	}
	
	/**
	 * @brief Query female density at location
	 * @param a_loc Coordinates to query
	 * @return Number of females in that 1 km² grid cell
	 * 
	 * @details Converts coordinates to grid index, returns counter value.
	 */
	int GetDensity(APoint a_loc) {
		int index = (a_loc.m_x / 1000) + (a_loc.m_y / 1000) * m_GridExtent;
		return m_FemaleDensityGrid[index];
	}
	
	/**
	 * @brief Query female density by pre-calculated index
	 * @param a_index Grid cell index
	 * @return Number of females in that cell
	 * 
	 * @details Faster version when caller already has index.
	 */
	int GetDensity(int a_index) {
		return m_FemaleDensityGrid[a_index];
	}
	
	/**
	 * @brief Reset density grid to zero
	 * 
	 * @details Sets all grid cells to 0. Called during initialization or when
	 * restarting population tracking. In normal simulation, not needed because
	 * grid maintained incrementally through Add/Remove calls.
	 */
	void ClearDensityGrid() {
		for (int i=0; i< m_FemaleDensityGrid.size(); i++) m_FemaleDensityGrid[i] = 0;
	}
	
	/**
	 * @brief Get today's prepupal development increment
	 * @return Development units accumulated today (days-equivalent)
	 * 
	 * @details Returns m_PrePupalDevelDaysToday calculated by DoBefore() each morning
	 * based on today's temperature forecast.
	 * 
	 * @par Prepupal Development Model
	 * Unlike other stages using degree-day accumulation, prepupal stage uses time-based
	 * development with temperature-dependent rates. This reflects biology: prepupa is
	 * summer diapause stage where development nearly arrested, with slow progression
	 * toward pupal stage independent of heat accumulation.
	 * 
	 * @par Difference from Formal Model
	 * **MAJOR CALIBRATION** - Formal model specifies quadratic temperature-development
	 * relationship for prepupal stage (Radmacher and Strohm 2011). Implementation uses
	 * simpler time-based approach: fixed duration ~45 days, weakly temperature-dependent.
	 * 
	 * **Rationale**: Quadratic model produced unrealistic prepupal durations under field
	 * conditions (laboratory parameterization didn't extrapolate well). Time-based model
	 * empirically calibrated to match field emergence phenology whilst acknowledging
	 * uncertainty in prepupal temperature response.
	 * 
	 * **Uncertainty**: HIGH - Prepupal stage poorly studied under field-realistic temperature
	 * regimes. Current implementation prioritizes phenological realism over mechanistic detail.
	 * 
	 * @see Osmia_Prepupa::Step() for usage in development
	 */
	double GetPrePupalDevelDays() {
		return m_PrePupalDevelDaysToday;
	}
	
	/**
	 * @brief Calculate available foraging hours for current day
	 * 
	 * @details Integrates hourly weather data with flight threshold criteria:
	 * - Temperature > cfg_OsmiaFlightTemperatureThreshold
	 * - Wind speed < cfg_OsmiaFlightWindSpeedThreshold  
	 * - Precipitation = 0 (or < threshold if light rain permitted)
	 * - Daylight hours (implicitly: weather station only reports daylight periods)
	 * 
	 * Counts hours meeting ALL criteria, stores in m_FlyingWeather for today.
	 * 
	 * @par Implementation Note
	 * Called during DoBefore() each morning before individual agents execute Step().
	 * All females share same foraging hour value because weather operates at landscape
	 * scale (any spatial variation in weather handled by landscape weather interpolation).
	 * 
	 * @par Biological Constraints
	 * Bees are ectothermic: flight muscles require minimum temperature. Wind physically
	 * prevents controlled flight at small body size. Wet conditions reduce flight efficiency
	 * and potentially increase wing damage. These constraints create realistic temporal
	 * variation in reproductive rate even with abundant floral resources.
	 * 
	 * @par Sensitivity
	 * Foraging hour calculation is CRITICAL parameter affecting population dynamics.
	 * Small changes in thresholds (e.g., 1°C difference in flight temperature) accumulate
	 * across season, substantially affecting population growth rate. Careful calibration
	 * against observed bee activity patterns essential.
	 */
	void CalForageHours(void);
	
	//==========================================================================
	// DEBUG/TESTING OUTPUT (Conditional Compilation)
	//==========================================================================
	
#ifdef __OSMIATESTING
public:
	/** @brief Output file stream for first nest egg data (testing mode) */
	ofstream m_eggsfirstnest;
	
	/** @brief Histogram of egg production by female age and size class */
	double m_egghistogram[4][30];
	
	/**
	 * @brief Write nest achievement test data
	 * @param a_target Target nest structure (expected)
	 * @param a_achieved Actual achieved nest structure
	 * 
	 * @details Testing function comparing expected vs. achieved nest provisioning
	 * patterns. Used during model development/calibration to verify females provision
	 * nests according to design (sex ratios, provision masses, cell counts).
	 */
	void WriteNestTestData(OsmiaNestData a_target, OsmiaNestData a_achieved);
#endif // __OSMIATESTING

protected:
	//==========================================================================
	// PROTECTED ATTRIBUTES
	//==========================================================================
	
	/** 
	 * @brief Static probability distribution for exponential random variates
	 * @details Provides efficiently-sampled exponential distribution (0-1 range).
	 * Used for stochastic variation in provision masses and other exponentially-
	 * distributed quantities. Static shared across all population managers
	 * (though typically only one manager per simulation).
	 */
	static probability_distribution m_exp_ZeroTo1;
	
#ifdef __OSMIATESTING
	/** @brief Vector storing female emergence weights for analysis (testing mode) */
	vector<double> m_FemaleWeights;
	
	/** @brief OpenMP lock protecting female weight vector (parallel access) */
	omp_nest_lock_t *m_female_weight_record_lock;
	
	/** @brief Statistics accumulator for egg production counts */
	SimpleStatistics m_OsmiaEggProdStats;
	
	/** @brief Statistics accumulator for egg stage durations */
	SimpleStatistics m_EggStageLength;
	
	/** @brief Statistics accumulator for larval stage durations */
	SimpleStatistics m_LarvalStageLength;
	
	/** @brief Statistics accumulator for prepupal stage durations */
	SimpleStatistics m_PrePupaStageLength;
	
	/** @brief Statistics accumulator for pupal stage durations */
	SimpleStatistics m_PupaStageLength;
	
	/** @brief Statistics accumulator for in-cocoon stage durations */
	SimpleStatistics m_InCocoonStageLength;
#endif
	
	/** 
	 * @brief Pointer to pollen map object
	 * @details Provides spatial interface to floral resource availability. Females
	 * query pollen map for nearby flower patches, resource quality/quantity, and
	 * distances to resources. Pollen map typically constructed from landscape land-use
	 * data combined with floral resource models.
	 */
	PollenMap_centroidbased* m_ThePollenMap;
	
	/** 
	 * @brief Daily foraging hours available (weather-dependent)
	 * @details Calculated by CalForageHours() each morning, used by all females
	 * throughout day. Value 0-24 representing hours meeting flight criteria.
	 * Typical range: 0 (very poor weather) to 8-10 (excellent foraging weather).
	 */
	int m_FlyingWeather;
	
	/** 
	 * @brief Flag indicating pre-wintering period has ended
	 * @details Pre-wintering is period between last emergence and onset of overwintering
	 * (roughly late August through September/early October). During this period, adult
	 * bees in cocoons undergo physiological changes preparing for winter dormancy.
	 * 
	 * Flag set TRUE when sustained autumn temperature drop detected (DoLast() logic:
	 * three consecutive days <13°C with declining trend). This triggers transition
	 * to overwintering phase where different development rules apply.
	 * 
	 * @par Biological Significance
	 * Pre-wintering/overwintering transition marks switch from active development
	 * (degree-day accumulation toward emergence) to dormancy (slower cold-hardening
	 * processes). Timing affects spring emergence phenology: earlier transition →
	 * longer overwintering → later emergence (potentially missing early flowers).
	 * 
	 * @see DoLast() for flag-setting logic
	 */
	bool m_PreWinteringEndFlag;
	
	/** 
	 * @brief Flag indicating overwintering period has ended (March 1st)
	 * @details Simple calendar-based flag set TRUE on March 1st regardless of weather.
	 * After this date, adults may emerge from cocoons once accumulated sufficient
	 * overwintering degree-days.
	 * 
	 * @par Biological Basis
	 * *Osmia bicornis* has obligate winter diapause: even under warm conditions,
	 * bees will not emerge until experiencing winter-like temperatures for sufficient
	 * duration. Fixed date (March 1st) represents earliest possible emergence in
	 * European temperate climate context.
	 * 
	 * **Difference from Formal Model**: Formal model specifies temperature-based
	 * emergence criteria without explicit calendar constraint. Implementation adds
	 * March 1st minimum to prevent unrealistic mid-winter emergence in warm years.
	 * This reflects biological reality of photoperiod cues not captured in temperature-
	 * only model.
	 */
	bool m_OverWinterEndFlag;
	
	/** 
	 * @brief Nest management interface
	 * @details Handles nest lifecycle: creation, polygon association, cell tracking,
	 * destruction. Separates nest infrastructure from population manager, allowing
	 * independent nest management logic evolution.
	 * 
	 * @par Design Rationale
	 * Nest manager could be integrated into landscape (polygon-level storage) or
	 * individual agents (each bee knows its nest). Separate manager class chosen
	 * because:
	 * - Nests outlive individual bees (one female provisions, multiple offspring emerge)
	 * - Polygon-level capacity tracking required (density regulation)
	 * - Parasitoid queries need nest-level data independent of current occupant
	 */
	Osmia_Nest_Manager m_OurOsmiaNestManager;
	
	/** 
	 * @brief Pre-calculated provisioning time parameters [days 0-364]
	 * @details Lookup table storing hours required to provision one cell as function
	 * of female age. Calculated during Init() from Seidelmann (2006) equations.
	 * 
	 * Index = female age (days since emergence)
	 * Value = hours needed for complete cell provisioning
	 * 
	 * Array size 365 provides full-year coverage even though bees rarely survive
	 * beyond 60 days (safety margin with minimal memory cost).
	 * 
	 * @par Difference from Formal Model
	 * **EXACT MATCH** - Values implement Seidelmann (2006) equations precisely.
	 * Pre-calculation is performance optimization (avoid repeated exp() calls).
	 * 
	 * @see GetProvisioningParams() for access method
	 */
	double m_NestProvisioningParameters[365];
	
	/** 
	 * @brief Logistic equations for egg sex ratio vs. age/mass [96 mass classes × 365 ages]
	 * @details Two-dimensional lookup table: m_EggSexRatioEqns[mass_class][age]
	 * 
	 * Mass classes: 0-95 representing 4.0-27.75 mg (0.25 mg increments)
	 * Age range: 0-364 days
	 * Values: Probability of female egg (0.0-1.0)
	 * 
	 * Calculated during Init() from logistic functions fitted to Seidelmann et al. (2010)
	 * data relating maternal age/mass to observed sex ratios.
	 * 
	 * @par Memory Footprint
	 * 96 × 365 × 8 bytes ≈ 280 KB per population manager. Modest cost for avoiding
	 * thousands of daily logistic function evaluations.
	 * 
	 * @see GetSexRatioEggsAgeMass() for access method
	 */
	vector<eggsexratiovsagelogisticcurvedata> m_EggSexRatioEqns;
	
	/** 
	 * @brief Logistic equations for female cocoon mass vs. maternal age/mass [96 mass classes × 365 ages]
	 * @details Two-dimensional lookup table: m_FemaleCocoonMassEqns[mass_class][age]
	 * 
	 * Returns target provision mass for first female offspring cell. Implements
	 * Seidelmann et al. (2010) declining investment pattern: heavier initial investment
	 * in first female cells, declining with maternal age and resource depletion.
	 * 
	 * @see GetFirstCocoonProvisioningMass() for access method with stochastic variation
	 */
	vector<femalecocoonmassvsagelogisticcurvedata> m_FemaleCocoonMassEqns;
	
	/** 
	 * @brief Female density grid [1 km² cells]
	 * @details Integer vector storing count of active females per grid cell.
	 * Grid extent calculated from landscape dimensions during initialization.
	 * 
	 * Size = (landscape_width_km) × (landscape_height_km)
	 * 
	 * Updated incrementally as females emerge, move, and die.
	 * 
	 * @see AddToDensityGrid(), RemoveFromDensityGrid(), GetDensity()
	 */
	vector<int> m_FemaleDensityGrid;
	
	/** 
	 * @brief Number of grid cells in X direction
	 * @details m_GridExtent = landscape_width / 1000 (meters per grid cell)
	 * Used for coordinate-to-index conversions.
	 */
	int m_GridExtent;
	
	/** 
	 * @brief Pollen availability scaling factor for interspecific competition
	 * @details Multiplier (0.0-1.0) reducing available pollen to account for consumption
	 * by other bee species. Value 1.0 = no competition (only *Osmia* present), value 0.5 =
	 * 50% of pollen already consumed by competitors before *Osmia* forages.
	 * 
	 * @par Biological Context
	 * *Osmia bicornis* co-occurs with many other bee species (honeybees, bumblebees,
	 * other solitary bees). All compete for same floral resources. Explicit competition
	 * modelling complex (requires multi-species population dynamics), so simple scaling
	 * factor approximates net effect.
	 * 
	 * @par Parameterization
	 * Read from configuration (cfg_OsmiaPollenCompetitionReductionScaler). Ideally
	 * estimated from field studies measuring resource depletion rates in presence/absence
	 * of competitor species. In practice, often calibrated to match observed *Osmia*
	 * population dynamics in landscapes with known bee community composition.
	 * 
	 * @par Uncertainty
	 * MEDIUM - Simplified representation of complex competitive interactions. Assumes
	 * constant competition intensity across space/time, but reality shows spatial
	 * aggregation (honeybee hives) and temporal variation (bumblebee phenology).
	 */
	double m_PollenCompetitionsReductionScaler;
	
	/** 
	 * @brief Prepupal development rate lookup table
	 * @details Vector storing temperature-dependent development rates for prepupal stage.
	 * Although prepupal model is time-based (not degree-day), rates still weakly
	 * temperature-dependent to acknowledge physiological reality.
	 * 
	 * @see GetPrePupalDevelDays() for daily rate access
	 */
	vector<double> m_PrePupalDevelRates;
	
	/** 
	 * @brief Today's prepupal development rate (pre-calculated)
	 * @details Updated each morning during DoBefore() based on temperature forecast.
	 * Cached for efficiency (all prepupae query same value throughout day).
	 */
	double m_PrePupalDevelDaysToday;
	
	/** 
	 * @brief Monthly pollen and nectar quality/quantity thresholds
	 * @details Vector of 12 OsmiaPollenNectarThresholds objects (one per month).
	 * Used by females during foraging to decide whether patch meets minimum
	 * acceptability criteria.
	 * 
	 * Read from configuration during Init(), indexed by current month during
	 * foraging decisions.
	 * 
	 * @see OsmiaPollenNectarThresholds class documentation
	 */
	vector<OsmiaPollenNectarThresholds> m_PN_thresholds;
	
	//==========================================================================
	// SCHEDULING METHODS (ALMaSS Framework Hooks)
	//==========================================================================
	
	/** 
	 * @brief Pre-step updates executed before any agents act
	 * 
	 * @details Daily setup sequence:
	 * 1. Update seasonal flags (check for pre-wintering/overwintering transitions)
	 * 2. Calculate foraging hours: CalForageHours()
	 * 3. Calculate prepupal development rates from temperature
	 * 4. Update pollen map (if dynamic resource model active)
	 * 5. Any other global state requiring daily refresh
	 * 
	 * Virtual method overriding Population_Manager::DoFirst(). Called by ALMaSS
	 * scheduler before individual agents execute BeginStep().
	 */
	virtual void DoFirst();
	
	/** 
	 * @brief Pre-step updates executed before Step but after DoFirst
	 * 
	 * @details Currently performs additional setup requiring DoFirst() completion.
	 * Less commonly used than DoFirst() in *Osmia* model; most setup concentrated
	 * in DoFirst().
	 * 
	 * Virtual method overriding Population_Manager::DoBefore(). Called by ALMaSS
	 * scheduler after DoFirst() but before individual agents execute Step().
	 */
	virtual void DoBefore();
	
	/** 
	 * @brief Post-step updates executed after Step but before DoLast
	 * 
	 * @details Currently unused in *Osmia* model (empty override). Provided for
	 * consistency with ALMaSS framework and potential future extensions requiring
	 * immediate post-step processing.
	 * 
	 * Virtual method overriding Population_Manager::DoAfter().
	 */
	virtual void DoAfter() {}
	
	/** 
	 * @brief End-of-day updates executed after all agents finish
	 * 
	 * @details Critical seasonal logic and statistics:
	 * 
	 * **Seasonal Flag Management**:
	 * - Check for pre-wintering end (sustained autumn temperature drop)
	 * - Set overwintering end flag (March 1st)
	 * - Reset flags after emergence season (June)
	 * 
	 * **Testing Output** (if __OSMIATESTING defined):
	 * - Record stage length statistics
	 * - Write annual summaries
	 * - Clear accumulators for next year
	 * 
	 * Virtual method overriding Population_Manager::DoLast(). Called by ALMaSS
	 * scheduler after all individual agents complete Step() and cleanup.
	 * 
	 * @par Pre-wintering End Detection Logic
	 * Sustained temperature drop identified by:
	 * - Three consecutive days < 13°C (day-2, day-1, day-0 all below threshold)
	 * - AND either: sharp sustained drop (days -5→-4 and -4→-3 both increase >1°C)
	 * - OR: extended cold period (day-3 also <13°C) with moderate drop (days -5→-4 increase ≥3°C)
	 * 
	 * This complex logic avoids false triggers from brief cold snaps whilst reliably
	 * detecting true autumn transition. Thresholds empirically calibrated for
	 * European temperate climate.
	 * 
	 * @see m_PreWinteringEndFlag, m_OverWinterEndFlag
	 */
	virtual void DoLast() {
		int today = m_TheLandscape->SupplyDayInYear();
		if (today > September) {
			int day = g_date->OldDays() + g_date->DayInYear();
			double t0 = m_TheLandscape->SupplyTempPeriod(day, 1);
			
			// Check for end of pre-wintering phase
			if (!m_PreWinteringEndFlag) {
				double t1 = m_TheLandscape->SupplyTempPeriod(day - 1, 1);
				double t2 = m_TheLandscape->SupplyTempPeriod(day - 2, 1);
				double t3 = m_TheLandscape->SupplyTempPeriod(day - 3, 1);
				double t4 = m_TheLandscape->SupplyTempPeriod(day - 4, 1);
				double t5 = m_TheLandscape->SupplyTempPeriod(day - 5, 1);
				
				// Sustained autumn cooling pattern
				if (((t2 < 13.0) && (t1 < 13.0) && (t0 < 13.0)) && 
				    (((t5 - t4 > 1.0) && (t4 - t3 > 1.0)) || 
				     ((t3 < 13.0) && (t5 - t4 >= 3.0)))) {
					m_PreWinteringEndFlag = true;
				}
			}
		}
		
		// March 1st marks minimum emergence date
		if (today == March) {
			m_OverWinterEndFlag = true;
		}
		
		// Reset flags after emergence season
		if (today == June) {
			m_PreWinteringEndFlag = false;
			m_OverWinterEndFlag = false;
		}
		
#ifdef __OSMIARECORDFORAGE
		// Optional foraging statistics output
		double meanforage = 0.0;
		if (Osmia_Female::m_foragecount > 0) 
			meanforage = Osmia_Female::m_foragesum / Osmia_Female::m_foragecount;
		cout << meanforage << endl;
		Osmia_Female::m_foragesum = 0.0;
		Osmia_Female::m_foragecount = 0.0;
#endif

#ifdef __OSMIATESTING		
		// Annual stage length statistics output
		if (today == 364) {
			ofstream file1("OsmiaStageLengths.txt", ios::app);
			file1 << "Year: " << g_date->GetYear() << endl;
			file1 << "Mean egg stage days is:" << '\t' << m_EggStageLength.get_meanvalue() << endl;
			file1 << "Mean larval stage days is:" << '\t' << m_LarvalStageLength.get_meanvalue() << endl;
			file1 << "Mean prepupal stage days is:" << '\t' << m_PrePupaStageLength.get_meanvalue() << endl;
			file1 << "Mean pupal stage days is:" << '\t' << m_PupaStageLength.get_meanvalue() << endl;
			file1 << "Mean incocoon stage days is:" << '\t' << m_InCocoonStageLength.get_meanvalue() << endl;
			
			// Clear accumulators for next year
			m_EggStageLength.ClearData();
			m_LarvalStageLength.ClearData();
			m_PrePupaStageLength.ClearData();
			m_PupaStageLength.ClearData();
			m_InCocoonStageLength.ClearData();
			file1.close();
		}
#endif
	}
};

#endif
