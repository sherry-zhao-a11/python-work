# 🧛‍♂️ Vampire Hunting v1.4.4 — Logical Inference Framework

## 📘 Overview
**Vampire Hunting v1.4.4** is a logical inference system that reconstructs the spread of a simulated infection among a group of participants.  
Using text-based records of daily test results and evening contact groups, the program determines *when* infections occurred, *who* transmitted them, and *how* the contagion evolved over time.

The repository contains:
- **`contact.py`** — the complete Python implementation of the reasoning framework.  
- **`Vampire test.pdf`** — the coursework document defining the analytical problem, input structure, and evaluation criteria.

---

## 🎯 Core Objectives
- Parse raw text input into structured data describing participants, daily test results, and contact networks.  
- Construct temporal **vampire knowledge tables (`vks`)** encoding each participant’s state (`H`, `V`, `U`) over time.  
- Apply forward, backward, and overnight propagation rules to maintain logical consistency.  
- Identify **infection windows (`iw`)** where transformations could logically occur.  
- Determine **potential sires (`ps`)** responsible for transmitting vampirism.  
- Refine infection data through **cyclic analysis (`cyclic_analysis2`)**, producing a coherent, stable model of infection spread.

---

## 📂 File Descriptions
| File | Description |
|------|--------------|
| **`contact.py`** | Core program implementing all logical reasoning modules (Sections 2–21). Performs validation, parsing, propagation, and inference. |
| **`Vampire test.pdf`** | Coursework description and project specification outlining the scenario, logic, and requirements. |

---

## 🧩 Section Structure (Summary)
The program is organized into logical sections, each performing a specific reasoning task:

| Section | Description |
|:--|:--|
| **2** | Validate that the input file exists before parsing. |
| **3** | Parse participants, days, AM test results, and PM contact groups into structured Python objects. |
| **4** | Display parsed data for verification and readability. |
| **5** | Determine which contact group a participant belongs to at a given time (`contacts_by_time()`). |
| **6–7** | Initialize vampire knowledge tables (`vks`) and apply test results. |
| **8–10** | Propagate confirmed states (forward for vampires, backward for humans, and overnight). |
| **11** | Integrate contact-group reasoning to constrain infection possibilities. |
| **12–13** | Identify **infection windows (`iw`)** and determine **potential sires (`ps`)**. |
| **14–15** | Refine sire data and infection windows for better temporal precision. |
| **16** | Update global knowledge tables based on refined infection information. |
| **17** | Perform the first stage of cyclic analysis to reach logical stability. |
| **18–19** | Classify vampires into **originals**, **newborns**, and **vamps_unclear**; compute **sire sets (`ss`)**. |
| **20–21** | Detect hidden vampires and run full-system cyclic analysis to ensure final consistency. |

---

## ⚙️ How to Run
1. Prepare an input text file formatted according to **Vampire test.pdf** (containing participants, daily tests, and contact data).  
2. Place the input file in the same directory as `contact.py`.  
3. Run the program in the terminal:
   ```bash
   python contact.py <input_filename>
