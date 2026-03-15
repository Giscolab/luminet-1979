# Luminet‑1979 — Visualisation WebGL2 d’un disque d’accrétion

Application front-end (sans build npm) qui affiche en temps réel un disque d’accrétion autour d’un trou noir stylisé.

## Ce que fait le projet

- rendu **WebGL2** full-screen (`canvas#c`) avec shader fragment principal ;
- mode géodésique « Kerr-lite / Schwarzschild » avec intégration visuelle orientée temps réel ;
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
