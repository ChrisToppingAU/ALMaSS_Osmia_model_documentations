---
name: midox
description: Comprehensive guide for creating MIDox (Model Implementation Documentation with Doxygen) for ecological simulation models. Use when documenting ALMaSS models or other ecological simulation implementations following the MIDox standard for FESMJ publication. Provides complete workflow from code enhancement through to dual-repository deployment.
license: CC-BY-4.0
---

# MIDox Documentation Generator

Comprehensive protocol for creating publication-ready Model Implementation Documentation with Doxygen (MIDox) for ecological simulation models, particularly within the ALMaSS framework.

## When to Use This Skill

Use this skill when:
- Documenting ecological simulation model implementations for FESMJ publication
- Creating MIDox documentation for ALMaSS spider models or similar agent-based models
- Enhancing existing minimally documented C++ ecological models
- Preparing implementation documentation following the three-paper sequence (Formal Model Ã¢â€ â€™ MIDox Ã¢â€ â€™ Testing & Calibration)

## Core MIDox Principles

**Documentation standards:**
- UK English spelling throughout
- No em dashes (use hyphens or restructure sentences)
- Literature citations for all substantive biological claims
- Doxygen-compatible comment formatting
- Dual-output approach: GitHub Pages (interactive) + Zenodo (archival DOI)
- Species names in italics every occurrence (e.g., *Osmia bicornis*, *O. bicornis*)
- Author placeholders in initial drafts (to be filled manually with correct affiliations)
- Seamless section flow (no "Part X" or "Continued from" markers in final document)

**Biological interpretation:**
- Transform minimal code comments into comprehensive biological documentation
- Every parameter requires biological rationale, empirical sources, uncertainty assessment
- Implementation decisions justified with reference to ecological theory
- Explicit acknowledgement of simplifications and assumptions

**Three-paper integration:**
- MIDox builds on the Formal Model (conceptual foundation)
- MIDox feeds into Testing & Calibration (empirical validation)
- Avoid redundancy whilst maintaining coherence across papers

## MIDox Workflow

### Phase 1: Preparation (1-2 hours)

1. **Gather source materials:**
   - Original C++ source files (.h and .cpp)
   - Formal Model paper (conceptual foundation) - **Extract exact citation** (authors, year, DOI)
   - Key literature (empirical parameters, biological basis) - **Create citation list**
   - Existing documentation (if any)

2. **Extract and verify key citations BEFORE generating narrative:**
   
   Create a reference file with exact citations for:
   - **Formal Model:** Copy complete citation from PDF title page
   - **Major empirical sources:** Extract from original papers (Radmacher & Strohm, Sgolastra et al., etc.)
   - **Framework papers:** ALMaSS, ODD protocol references
   - **Key biology papers:** Species reviews, behavior studies
   
   **Format as a prompt-ready citation block:**
   ```
   VERIFIED CITATIONS TO USE:
   
   Formal Model:
   Ziółkowska E, Bednarska AJ, Laskowski R, Topping CJ (2023). The Formal Model 
   for the solitary bee Osmia bicornis L. agent-based model. Food and Ecological 
   Systems Modelling Journal 4: e102102. https://doi.org/10.3897/fmj.4.102102
   
   Key Empirical Sources:
   Radmacher S, Strohm E (2011). Effects of constant and fluctuating temperatures 
   on the development of the solitary bee Osmia bicornis (Hymenoptera: Megachilidae). 
   Apidologie 42(6): 711-720. https://doi.org/10.1007/s13592-011-0078-9
   
   [Add all other key papers...]
   ```
   
   Include this citation block in every prompt to AI to ensure accurate references.

2. **Review existing code:**
   - Identify all classes, methods, parameters
   - Map inheritance hierarchies
   - Note areas with minimal documentation
   - Check for `_enh` (enhanced) versions if updating

3. **Read MIDox standards:**
   - See [references/midox_appendices.md](references/midox_appendices.md) for complete section structure
   - Review example documentation patterns
   - Understand Doxygen hyperlinking conventions

### Phase 2: File-by-File Enhancement (30-60 hours for typical model)

**Enhancement strategy:**
- Work file-by-file due to token limits
- Start with base classes (most reusable across species)
- Follow with species-specific implementations
- End with population managers and utility functions

**For each file, document:**

1. **File header:**
   ```cpp
   /**
    * @file Erigone.cpp
    * @brief Implementation of Erigone atra (black money spider) life cycle
    * @details Comprehensive implementation including development, reproduction, 
    * mortality, and dispersal processes calibrated from field studies. See 
    * Thorbek & Topping (2005) for model foundation and De Keer & Maelfait (1988) 
    * for species biology.
    * @author Original: Topping & Thorbek. Enhanced: Your Name
    * @date Original: 2004. Enhanced: 2025
    * @ingroup Spider_Models
    */
   ```

2. **Class documentation:**
   ```cpp
   /**
    * @class Erigone
    * @brief Agent-based implementation of Erigone atra life cycle
    * @details This class implements individual spider lifecycle including:
    * - Temperature-dependent development using degree-day model
    * - Stage-specific mortality (eggs, juveniles, adults)
    * - Reproductive biology with pre-oviposition and egg production
    * - Habitat-dependent dispersal including ballooning
    * 
    * @par Biological Foundation
    * Implementation based on field studies by Thorbek (2002) conducted in 
    * Danish agricultural landscapes. Development parameters from De Keer & 
    * Maelfait (1988), dispersal from field observations showing strong habitat 
    * preferences for winter cereals.
    * 
    * @par Implementation Notes
    * Uses Sharpe-Schoolfield temperature response for development to handle
    * both cold and heat stress. Mortality implemented as daily probability 
    * rather than age-specific schedules for computational efficiency.
    */
   ```

3. **Method documentation:**
   ```cpp
   /**
    * @brief Calculate daily development progress for juvenile spiders
    * @details Implements Sharpe-Schoolfield temperature-dependent development 
    * model with low and high temperature inhibition. Development accumulates 
    * as fraction of total requirement (0.0 to 1.0).
    * 
    * @param temperature Current ambient temperature (Ã‚Â°C)
    * @return Daily development increment (fraction of total development)
    * 
    * @par Biological Rationale
    * Spider development ceases below threshold temperature (0Ã‚Â°C for E. atra)
    * and slows above optimal range due to enzyme denaturation. Field 
    * observations show juveniles require approximately 1800 degree-days 
    * for maturation (De Keer & Maelfait 1988).
    * 
    * @par Implementation Details
    * Uses pre-calculated rate at 25Ã‚Â°C (ERIGONE_JUVDEVELRH) and activation
    * energy (ERIGONE_JUVDEVELHA) to calculate temperature-specific rates.
    * Avoids recalculating exponentials on every call for efficiency.
    * 
    * @par Uncertainty
    * Laboratory-derived parameters may not fully capture field conditions.
    * Development threshold shows population variation (Ã‚Â±2Ã‚Â°C). Calibration 
    * adjusted DD requirement to 1800 from laboratory 1600 to match field 
    * emergence timing.
    * 
    * @see CalculateEggDevelopment() for similar implementation in egg stage
    * @see Thorbek & Topping (2005) for model calibration approach
    */
   double CalculateJuvenileDevelopment(double temperature);
   ```

4. **Parameter documentation:**
   ```cpp
   /**
    * @var ERIGONE_JUVDEVELCONSTTWO
    * @brief Total degree-days required for juvenile development to adult
    * @details Default: 1800 degree-days above 0Ã‚Â°C threshold
    * 
    * @par Empirical Basis
    * Based on De Keer & Maelfait (1988) laboratory study showing 1600 DD at 
    * 0Ã‚Â°C threshold. Value increased to 1800 DD during calibration to match 
    * field emergence timing observed by Thorbek (2002) in Danish winter wheat.
    * 
    * @par Biological Interpretation
    * Represents combined developmental requirements across 6-8 instars from 
    * hatching to sexual maturity. Higher field value likely reflects energetic 
    * costs of foraging and predator avoidance not present in laboratory.
    * 
    * @par Sensitivity
    * HIGH - Directly affects maturation timing and population phenology. Ã‚Â±10% 
    * change shifts peak adult abundance by 7-10 days in simulations.
    * 
    * @par Valid Range
    * [1500, 2200] degree-days. Values below 1500 produce unrealistically early 
    * emergence; above 2200 delays reproduction beyond observed field patterns.
    */
   cfg_erigone_juvdevelconsttwo ERIGONE_JUVDEVELCONSTTWO("ERIGONE_JUVDEVELCONSTTWO", CFG_CUSTOM, 1800);
   ```

**Quality indicators for enhanced files:**
- Every public method has biological rationale
- Every parameter has empirical source and uncertainty assessment
- Inheritance relationships clearly explained
- Complex algorithms have step-by-step commentary
- Simplifications and assumptions explicitly stated
- File typically 50-120% longer than original

### Phase 3: Narrative Documentation (8-12 hours)

Create comprehensive MIDox markdown document following structure in [references/midox_structure.md](references/midox_structure.md):

1. **Introduction** (~500 words): Position within three-paper sequence
2. **Model Overview** (~1500 words): Architecture, emergence, spatial/temporal scales
3. **Scheduling and Workflow** (~1000 words): Execution sequence, pseudocode
4. **Implementation Details** (~2000 words): Component descriptions with code links
5. **Parameters** (~1000 words): Organisation, sources, uncertainty
6. **Inputs** (~800 words): Required data, formats, preprocessing
7. **Outputs** (~800 words): Generated data, formats, interpretation
8. **Implementation Discussion** (~1500 words): Trade-offs, limitations, future directions

**Table formatting for Doxygen markdown:**
- All rows MUST have the same number of columns as the header
- Section headers within tables (e.g., "**Egg Stage**") must span all columns with empty cells
- Incorrect: `| **Egg Stage** |` (1 column)
- Correct: `| **Egg Stage** |||||` (6 columns for a 6-column table)
- Tables must be continuous - don't break into separate tables unnecessarily

**Final section (replaces appendices):**
9. **Documentation Access**: Brief section (~200 words) directing readers to:
   - Interactive documentation: GitHub Pages URL
   - Archived version: Zenodo DOI with badge
   - Source code repository: GitHub URL
   - Compilation instructions: Link to README in repository

**Integration with Doxygen:**

CRITICAL: Doxygen treats markdown files (.md) differently from C++ comment blocks. In pure markdown files, you CANNOT use bare @ref tags—they appear as literal text. You must use markdown link syntax.

**For class names in bullet lists:**
```markdown
* [Osmia_Egg](@ref Osmia_Egg): From laying to hatching (temperature-dependent duration)
* [Osmia_Larva](@ref Osmia_Larva): Feeding and cocoon construction (resource-dependent growth)
* [Osmia_Prepupa](@ref Osmia_Prepupa): Summer diapause (time-based, weakly temperature-dependent)
```

**For method names (always include parentheses):**
```markdown
* [DoFirst()](@ref DoFirst()): Update global environmental conditions
* [DoBefore()](@ref DoBefore()): Pre-step calculations
* [Step()](@ref Step()): Individual agents execute behaviour
```

**For method names IN TABLES (must use fully qualified names):**
```markdown
| [Osmia_Base](@ref Osmia_Base) | Purpose | [GetTemp()](@ref Osmia_Base::GetTemp()), [SetTemp()](@ref Osmia_Base::SetTemp()) |
```
**CRITICAL:** In tables, methods MUST use `ClassName::MethodName()` format:
- Wrong: `[Step()](@ref Step())` - won't link in tables
- Right: `[Step()](@ref Osmia_Female::Step())` - creates proper link

**For class members:**
```markdown
* [Osmia_Base::m_Age](@ref Osmia_Base::m_Age): Agent age in days
```

**In flowing prose:**
```markdown
All agent classes inherit from [Osmia_Base](@ref Osmia_Base), which provides common attributes (age, mass, location) and methods ([CheckMortality()](@ref CheckMortality()), temperature access).
```

**Why this is necessary:**
- In C++ comment blocks (`/** ... */`), bare `@ref ClassName` works
- In pure markdown files (.md), you MUST use `[ClassName](@ref ClassName)` syntax
- The markdown link format `[text](url)` wraps the @ref command
- Without brackets, Doxygen treats @ref as literal text in markdown

**Table formatting rules:**
- **CRITICAL:** All table rows must have identical column counts
- Section headers must span all columns: `| **Section Name** ||||...||`
- Count pipes in header row, match exactly in all data rows
- Wrong: `| **Egg Stage** |` (breaks table rendering)
- Right: `| **Egg Stage** |||||` (for 6-column table)

**General principles:**
- ALWAYS use markdown link syntax `[ClassName](@ref ClassName)` in .md files
- ALWAYS include parentheses for methods: `[MethodName()](@ref MethodName())`
- ALWAYS match column counts in table rows
- Test all links in generated HTML before publication

**TOKEN-MANAGED WORKFLOW FOR AI ASSISTANTS:**

When working with AI to create narrative documentation, use systematic part-by-part approach with checkpoint saves:

**Step 1: Create Part 1 (Sections 1-2) - Save immediately**
- Title block with PLACEHOLDER for authors (to be filled manually)
- Abstract with keywords
- Section 1: Introduction (full detail)
- Section 2: Model Overview (full detail)
- SAVE to outputs as Part1.md

**Step 2: Create Part 2 (Sections 3-4) - Save immediately**  
- Section 3: Scheduling and Workflow (full detail with pseudocode)
- Section 4: Implementation Details (full detail with algorithms)
- SAVE to outputs as Part2.md

**Step 3: Create Part 3 (Sections 5-6) - Save immediately**
- Section 5: Parameters (complete tables with all justifications)
- Section 6: Inputs (all data requirements)
- SAVE to outputs as Part3.md

**Step 4: Create Part 4 (Sections 7-9) - Save immediately**
- Section 7: Outputs (all output specifications)
- Section 8: Implementation Discussion (trade-offs, limitations)
- Section 9: Documentation Access (links to GitHub Pages/Zenodo)
- Complete References section
- SAVE to outputs as Part4.md

**Step 5: Combine all parts into single seamless document**
- Concatenate Part1.md + Part2.md + Part3.md + Part4.md
- REMOVE part heading markers (e.g., "# Part 2:", "Continued from Part X")
- Ensure smooth transitions between sections
- Apply consistent formatting throughout
- SAVE final combined document as ModelName_MIDox_Complete.md

**Step 6: Apply formatting requirements**
- Replace all instances of species names with italic markdown (e.g., "Osmia bicornis" â†’ "*Osmia bicornis*", "O. bicornis" â†’ "*O. bicornis*")
- Verify UK English spelling throughout
- Check no em dashes present
- Ensure all citations formatted consistently
- SAVE corrected version

**Token Budget Management:**
- Each part: ~20-30K tokens
- Total for complete narrative: ~80-120K tokens  
- Leaves buffer for revisions and formatting
- If hitting limits: save more frequently, create smaller parts

**Quality Checks Between Parts:**
- Each part is complete standalone section
- No "to be continued" or "TBD" placeholders
- All algorithms fully specified in pseudocode
- All parameters have values, sources, rationales
- Biological justifications complete

This approach ensures complete documentation even with token constraints, with each checkpoint providing usable content.

### Phase 4: Doxygen Configuration (2-3 hours)

Create or update Doxyfile:

```
PROJECT_NAME           = "Osmia bicornis Population Model"
PROJECT_BRIEF          = "Agent-based bee model within ALMaSS framework"
OUTPUT_DIRECTORY       = ./docs
INPUT                  = ./src ./docs/Osmia_MIDox_Complete.md
FILE_PATTERNS          = *.cpp *.h *.md
RECURSIVE              = YES
USE_MDFILE_AS_MAINPAGE = ./docs/Osmia_MIDox_Complete.md
GENERATE_HTML          = YES
GENERATE_LATEX         = NO
HTML_OUTPUT            = html
EXTRACT_ALL            = YES
EXTRACT_PRIVATE        = NO
EXTRACT_STATIC         = YES
SOURCE_BROWSER         = YES
INLINE_SOURCES         = NO
REFERENCED_BY_RELATION = YES
REFERENCES_RELATION    = YES
MARKDOWN_SUPPORT       = YES
AUTOLINK_SUPPORT       = YES
HAVE_DOT               = YES
CLASS_DIAGRAMS         = YES
CALL_GRAPH             = YES
CALLER_GRAPH           = YES
```

**CRITICAL for markdown integration:**
- `INPUT` must include both source directory AND the markdown file path
- `FILE_PATTERNS` must include `*.md` to process markdown files
- `USE_MDFILE_AS_MAINPAGE` makes the MIDox document the documentation homepage
- `MARKDOWN_SUPPORT = YES` enables markdown processing (should be default but be explicit)
- `AUTOLINK_SUPPORT = YES` enables @ref tag processing

Test generation locally:
```bash
doxygen Doxyfile
cd docs/html
python3 -m http.server 8000
# View at localhost:8000
```

### Phase 5: Repository Setup and Deployment (3-5 hours)

**GitHub repository structure:**
```
erigone-model/
Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ src/                          # Enhanced source code
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Erigone.h
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Erigone.cpp
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Erigone_Population_Manager.h
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Erigone_Population_Manager.cpp
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Spider_BaseClasses.h
Ã¢â€â€š   Ã¢â€â€Ã¢â€â‚¬Ã¢â€â‚¬ Spider_BaseClasses.cpp
Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ docs/
Ã¢â€â€š   Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Erigone_MIDox.md         # Narrative documentation
Ã¢â€â€š   Ã¢â€â€Ã¢â€â‚¬Ã¢â€â‚¬ html/                     # Generated by Doxygen
Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ Doxyfile                      # Doxygen configuration
Ã¢â€Å“Ã¢â€â‚¬Ã¢â€â‚¬ README.md                     # Overview with deployment links
Ã¢â€â€Ã¢â€â‚¬Ã¢â€â‚¬ LICENSE                       # Open source licence

```

**Deployment steps:**

1. **Generate HTML:**
   ```bash
   doxygen Doxyfile
   ```

2. **Deploy to GitHub Pages:**
   - Push to GitHub repository
   - Enable GitHub Pages in repository settings
   - Select 'docs/html' as source directory
   - Note URL: `https://username.github.io/erigone-model/`

3. **Archive to Zenodo:**
   - Create new Zenodo deposit
   - Upload complete repository as .zip
   - Add comprehensive metadata (authors, description, keywords)
   - Publish to receive DOI
   - Note DOI: `10.5281/zenodo.XXXXXXX`

4. **Update README with links and add final Documentation Access section to MIDox:**

   In MIDox narrative final section:
   ```markdown
   ## 9. Documentation Access
   
   Complete interactive documentation with full code reference and API details is available online:
   
   **Interactive documentation:** [https://username.github.io/model-name/](https://username.github.io/model-name/)
   
   **Archived version with DOI:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   
   **Source code repository:** [https://github.com/username/model-name](https://github.com/username/model-name)
   
   The GitHub Pages site provides:
   - Complete API documentation generated by Doxygen
   - Searchable class and method references  
   - Inheritance diagrams and call graphs
   - Cross-referenced source code
   - Detailed parameter documentation with @ref links to code
   
   The Zenodo archive provides a permanent, citable snapshot with DOI for long-term reproducibility.
   
   For compilation instructions, system requirements, and configuration examples, see README.md in the source repository.
   ```

   In repository README.md:
   ```markdown
   # Model Name - MIDox Documentation
   
   ## Access Documentation
   
   **Interactive documentation:** https://username.github.io/model-name/
   
   **Archived version:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   
   ## Citation
   
   If using this model, please cite:
   
   [Author names]. ([Year]). [Model name]: MIDox implementation documentation. 
   Food and Ecological Systems Modelling Journal, XX, XXX-XXX. DOI: 10.XXXX/fesmj.XXXX
   ```

### Phase 6: Quality Assurance (2-4 hours)

**Documentation review checklist:**
- [ ] 100% function documentation coverage
- [ ] All parameters have biological rationale
- [ ] All parameters cite empirical sources
- [ ] **CRITICAL: Every citation verified against original source**
- [ ] **CRITICAL: Formal model citation matches PDF exactly**
- [ ] **CRITICAL: All DOIs tested and resolve correctly**
- [ ] Uncertainty explicitly assessed for key parameters
- [ ] Inheritance relationships clearly explained
- [ ] Complex algorithms have biological interpretation
- [ ] Doxygen generates without errors/warnings
- [ ] All hyperlinks functional (narrative Ã¢â€ â€™ code Ã¢â€ â€™ narrative)
- [ ] **CRITICAL: All tables render correctly (check column counts match)**
- [ ] README includes both GitHub Pages and Zenodo links
- [ ] Narrative document follows MIDox section structure
- [ ] UK English spelling throughout
- [ ] No em dashes in generated text
- [ ] File names consistent with conventions

**Independent review:**
Have someone unfamiliar with the model:
- Navigate the documentation
- Attempt to understand a specific feature
- Verify they can trace from narrative to code and back
- Check for broken links or unclear explanations

## Time Investment Estimates

For a model of similar complexity to Erigone (7 C++ files, ~3500 lines code):

**File enhancement:** 30-60 hours
- Base classes: 8-12 hours
- Species-specific: 12-20 hours
- Population manager: 6-10 hours
- Utility functions: 4-8 hours

**Narrative documentation:** 8-12 hours

**Technical setup:** 5-8 hours

**Review and refinement:** 4-7 hours

**Total:** 47-87 hours

Subsequent models reusing base classes reduce time by ~30%.

## Extension to New Models

When documenting additional species (e.g., Oedothorax after Erigone):

1. **Leverage base class documentation** - Already complete, reusable
2. **Focus on species-specific differences** - What makes this species unique?
3. **Copy-adapt documentation patterns** - Use Erigone as template
4. **Update shared components minimally** - Only if genuinely improved
5. **Maintain consistency** - Same terminology, citation style, structure

Time investment typically 60-70% of first implementation.

## Common Pitfalls and Solutions

**Problem:** Documentation becomes disconnected from code
**Solution:** Use Doxygen's `@ref` tags liberally; test all links before publication

**Problem:** Biological rationale too vague
**Solution:** Ask "Why this value?" and "What if it were different?" for every parameter

**Problem:** Uncertainty assessment superficial
**Solution:** Explicitly state: measurement error, inter-population variation, lab vs field differences

**Problem:** Enhanced files not synchronised with project directory
**Solution:** Use `_enh` suffix; clearly mark which are working versions

**Problem:** Token limits prevent complete file enhancement
**Solution:** Work file-by-file; save progress; return to complete next file

**Problem:** Doxygen generation fails
**Solution:** Check for unmatched comment delimiters (`/**` without `*/`); validate with `doxygen -g test.cfg`

## AI Interaction Tips

When working with AI assistants on MIDox generation:

1. **Provide complete context upfront:**
   - Upload MIDox paper and appendices
   - Share Formal Model paper
   - Include key biological literature
   - Show enhanced file examples

2. **Work incrementally:**
   - Complete one file fully before next
   - Review AI-generated documentation carefully
   - Verify biological accuracy (AI may confabulate sources)

3. **Be specific about requirements:**
   - "Add biological rationale for this parameter"
   - "Explain why we use daily mortality not age-specific"
   - "Cite the empirical source for this threshold"

4. **Iterative refinement:**
   - First pass: structure and coverage
   - Second pass: biological depth
   - Third pass: uncertainty and limitations

5. **Verify all citations - CRITICAL:**
   
   **AI assistants frequently hallucinate or construct incorrect citations.** This is one of the most common and serious errors in AI-generated scientific writing.
   
   **For EVERY citation in the generated narrative:**
   - [ ] Verify the paper actually exists (search DOI or title)
   - [ ] Check authors are correct and in correct order
   - [ ] Verify publication year is accurate
   - [ ] Confirm journal name, volume, and page numbers
   - [ ] Test that DOI links to the correct paper
   - [ ] Verify claimed findings match what the paper actually says
   
   **Common AI citation errors:**
   - Wrong year (inferring from context rather than checking)
   - Incorrect author lists (confusing similar papers)
   - Hallucinated DOIs that look plausible but don't exist
   - Mixing authors from different papers
   - Claiming findings the paper doesn't actually contain
   
   **Safe citation workflow:**
   1. Provide AI with exact citations from source PDFs (copy-paste)
   2. For key papers (formal model, major empirical sources), paste full citation in prompt
   3. After generation, verify every single citation against original sources
   4. Never trust AI-constructed citations without verification
   
   **Special attention needed for:**
   - The formal model citation (see References section below)
   - Empirical parameter sources (Radmacher & Strohm, Sgolastra et al., etc.)
   - Methodology papers (ALMaSS framework, ODD protocol)
   - Species biology reviews

## Publication Integration

**FESMJ submission:**
- MIDox narrative as main manuscript body
- Parameter tables as supplementary material
- GitHub Pages URL in Methods section
- Zenodo DOI in Data Availability statement

**Long-term maintenance:**
- GitHub Pages documentation is "living" (updates with code)
- Zenodo archive is permanent (snapshot with DOI)
- Major updates create new Zenodo versions with new DOIs
- Original MIDox paper remains unchanged; readers directed to current docs

## References for MIDox Implementation

**CRITICAL - Citation Accuracy for ALL References:**

AI assistants frequently generate plausible-looking but completely incorrect citations. This is one of the most serious and common errors in AI-generated scientific documentation. **Never trust AI-constructed citations without verification.**

**Mandatory verification for every citation:**
1. Extract EXACT citations from source PDFs (title pages, not inferred)
2. Verify author list, year, journal volume, page numbers, and DOI match source exactly
3. Test that DOI actually resolves to the correct paper
4. Confirm claimed findings match what the paper actually says
5. Do NOT allow AI to construct citations from memory or inference

**Common citation errors to watch for:**
- Wrong publication year (AI infers from context or confuses with other papers)
- Incorrect author lists (mixing authors from different papers, wrong order)
- Hallucinated DOIs that look plausible but don't exist or link to wrong papers
- Wrong journal volumes or page numbers
- Misattributed findings (claiming Paper A's results are from Paper B)

**Required workflow:**
1. **Before generating narrative:** Create verified citation list from source PDFs
2. **Include in prompts:** Paste exact citations into every AI prompt
3. **After generation:** Verify EVERY citation against original sources
4. **Fix immediately:** Correct any errors before finalizing document

**Examples of correct citation formats:**

Formal Model:
```
Ziółkowska E, Bednarska AJ, Laskowski R, Topping CJ (2023). The Formal Model for 
the solitary bee Osmia bicornis L. agent-based model. Food and Ecological Systems 
Modelling Journal 4: e102102. https://doi.org/10.3897/fmj.4.102102
```

Empirical Study:
```
Radmacher S, Strohm E (2011). Effects of constant and fluctuating temperatures on 
the development of the solitary bee Osmia bicornis (Hymenoptera: Megachilidae). 
Apidologie 42(6): 711-720. https://doi.org/10.1007/s13592-011-0078-9
```

Framework Paper:
```
Topping CJ, Høye TT, Olesen CR (2003). Opening the black box—Development, testing 
and documentation of a mechanistically rich agent-based model. Ecological Modelling 
166(1-2): 109-123. https://doi.org/10.1016/S0304-3800(03)00140-5
```

**Why this matters:**
Incorrect citations undermine scientific credibility, prevent readers from finding sources, attribute findings incorrectly, and can constitute scientific misconduct if errors appear systematic.

See [references/midox_literature.md](references/midox_literature.md) for:
- Complete MIDox editorial paper
- MIDox appendices with section structure
- Example implementations (hare model)
- Doxygen manual and best practices
- ALMaSS framework documentation

## Support Files

This skill includes:

- **references/midox_appendices.md** - Complete section-by-section structure
- **references/midox_structure.md** - Detailed content requirements
- **references/doxygen_patterns.md** - Comment formatting examples
- **references/midox_literature.md** - Key papers and resources
- **examples/enhanced_files/** - Example enhanced C++ files
- **templates/Doxyfile** - Standard Doxygen configuration
- **templates/README_template.md** - Repository README structure

Read these as needed using the `view` tool.

**Citation verification (mandatory):**
For each citation in the narrative:
- [ ] Authors match source exactly (order, spelling, initials)
- [ ] Publication year is correct
- [ ] Journal name, volume, and page numbers are accurate
- [ ] DOI resolves to the correct paper
- [ ] Quoted or paraphrased findings actually exist in the paper
- [ ] No citations are AI hallucinations or confabulations
