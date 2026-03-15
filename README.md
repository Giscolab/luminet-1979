# Luminet‑1979 — Visualisation WebGL2 d’un disque d’accrétion

Application front-end (sans build npm) qui affiche en temps réel un disque d’accrétion autour d’un trou noir stylisé.

## Ce que fait le projet

- rendu **WebGL2** full-screen (`canvas#c`) avec shader fragment principal ;
- mode géodésique « Kerr-lite / Schwarzschild » avec intégration visuelle orientée temps réel ;
- nouveau mode de physique des rayons (`effective` rapide vs `kerr_full` basé sur constantes conservées Kerr `E/Lz/Q`) piloté par un toggle UI ;
- contrôles interactifs (inclinaison, rayon externe, vitesse, spin, lensing, qualité, exposition, etc.) ;
- caméra orbitale ou manuelle ;
- modes d’affichage `NORMAL` et `BLOOM` ;
- émissivité disque `NT` (par défaut) ou `legacy stylized`.

## Nouveautés (mise à jour)

- **Persistance locale des réglages** : les paramètres UI sont sauvegardés dans `localStorage` puis restaurés au rechargement.
- **Raccourci `R`** : réinitialise les paramètres aux valeurs par défaut puis resauvegarde l’état.
- Le README reflète désormais l’état réel de `script.js` (il n’est plus « non modifié »).

## Fichiers

- `index.html` : structure de la page, overlays et contrôles
- `style.css` : apparence des panneaux/contrôles
- `script.js` : moteur WebGL2, shaders, interactions utilisateur, persistance locale

## Lancer le projet

Aucune dépendance à installer.

### Option 1 — ouverture directe
Ouvrir `index.html` dans un navigateur compatible WebGL2.

### Option 2 — serveur statique local
```bash
python3 -m http.server 8000
```
Puis ouvrir : <http://localhost:8000/>

## Contrôles clavier

- `A` : badge analytique
- `G` : badge géodésique RK4
- `L` : badge LUT expérimental
- `R` : reset des paramètres (avec sauvegarde)

## Remarques

- Projet orienté rendu visuel temps réel (pas solveur GR complet).
- Un GPU/driver récent améliore nettement la fluidité.


## Compromis performance / précision (physicsMode)

- `effective` (défaut) : conserve l'accélération géométrique historique (`geoAccel`), très stable et rapide pour l'aperçu temps réel.
- `kerr_full` : active une voie RK4 dédiée exploitant une formulation Kerr en coordonnées de type Boyer-Lindquist avec invariants (`E`, `Lz`, `Q`) pour mieux respecter la géométrie proche horizon.
- Impact attendu : `kerr_full` augmente le coût par pas (plus de calculs transcendants + potentiels radiaux/polaires), donc fps plus faible surtout en qualité élevée.
- Garde-fous numériques ajoutés : pas minimal adaptatif, clamp robuste autour de l'horizon et seuil de sortie en cas d'état numérique divergent.
- Recommandation pratique : utiliser `effective` pour l'exploration interactive, puis basculer en `kerr_full` pour inspecter les trajectoires critiques / proches horizon.
