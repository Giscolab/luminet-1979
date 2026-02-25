# Disque d’accrétion relativiste — Visualisation WebGL2

Cette application est une **visualisation temps réel** d’un disque d’accrétion autour d’un trou noir, avec un mode géodésique RK4, des réglages caméra/physique et une interface d’instrumentation minimaliste.

> ⚠️ Contrainte respectée : le fichier `script.js` n’a pas été modifié.

---

## 1) Analyse complète du projet (état actuel)

## Objectif fonctionnel
Le projet propose une scène interactive WebGL2 qui simule :
- la courbure des rayons lumineux,
- les variations de luminosité/couleur (effets Doppler/relativistes simplifiés),
- un disque d’accrétion paramétrable,
- plusieurs modes d’affichage (normal / bloom, badge A/G/L).

## Architecture (simple et lisible)
Le dépôt contient 3 fichiers applicatifs + ce README :
- `index.html` : structure de l’interface (canvas + overlays + contrôles).
- `style.css` : habillage visuel et ergonomie des panneaux de contrôle.
- `script.js` : moteur WebGL2 (shaders, boucle de rendu, interactions).

## Détail de l’interface HTML
La page contient :
- un `<canvas id="c">` qui porte le rendu principal ;
- une barre de progression LUT (`#lut-progress`) ;
- un réticule dynamique (`#crosshair`) ;
- un rappel de raccourcis clavier (`.hint`) ;
- trois panneaux UI :
  - statistiques (inclinaison, rayon externe, fps),
  - contrôles sliders/checkbox,
  - boutons de mode (`NORMAL` / `BLOOM`) ;
- un badge de mode (`#mode-badge`).

## Moteur JavaScript (non modifié)
Le script :
- initialise WebGL2 et compile les shaders ;
- envoie des uniforms de simulation (`iIncl`, `iRout`, `iSpin`, `iQuality`, etc.) ;
- trace les rayons avec intégration RK4 et logique de disque/horizon ;
- applique colorimétrie, glow/fog et bloom optionnel ;
- gère les interactions utilisateur (sliders, touches `A/G/L`, mode bloom, caméra orbitale/libre, crosshair, FPS).

## Contrôles disponibles
- **Clavier** : `A`, `G`, `L` pour les badges/modes signalés.
- **Souris / tactile** : interaction caméra / rayon externe selon logique JS.
- **Sliders** : inclinaison, rayon externe, vitesse rotation, yaw/pitch/dist caméra, épaisseur disque, lensing, spin, qualité.
- **Toggle** : suivi d’orbite caméra (active/désactive sliders manuels caméra).
- **Boutons** : `NORMAL` et `BLOOM`.

## État qualité observé avant refonte CSS
- Interface déjà fonctionnelle mais visuellement hétérogène sur certains panneaux.
- Contraste parfois limite sur labels secondaires.
- Responsive incomplet pour petits écrans (densité de contrôles élevée).
- États UI (hover, active, disabled) perfectibles pour le confort visuel.

---

## 2) Refonte CSS réalisée (sans toucher au JavaScript)

Le fichier `style.css` a été modernisé pour améliorer le confort visuel et la cohérence avec l’UI existante :

- **Système de variables CSS (`:root`)** pour homogénéiser palette, contrastes, rayons, ombres.
- **Panels glassmorphism léger** (fond translucide + blur + bordures douces) pour lisibilité sans masquer la scène.
- **Overlay en grille responsive** pour une meilleure adaptation desktop/tablette/mobile.
- **Typographies stabilisées** (textes UI + monospace pour valeurs numériques).
- **Sliders améliorés** : piste, thumb, états hover + état disabled (utile quand caméra orbitale est active).
- **Boutons de mode renforcés** : états hover/active plus lisibles.
- **Mode badge et hint** harmonisés (contraste et lisibilité accrue).
- **Crosshair** plus net avec halo subtil.
- **Media queries** dédiées (`980px`, `600px`) pour éviter la surcharge visuelle sur petit écran.

Résultat : interface plus propre, plus lisible, plus cohérente avec le rendu astrophysique et les interactions déjà présentes.

---

## 3) Lancer le projet

Aucune dépendance npm :

1. Ouvrir `index.html` dans un navigateur compatible WebGL2.
2. Ou lancer un serveur local statique :

```bash
python3 -m http.server 8000
```

Puis ouvrir : `http://localhost:8000/`.

---

## 4) Compatibilité recommandée

- Navigateurs récents : Chrome / Edge / Firefox récents (WebGL2 requis).
- GPU/driver à jour conseillé pour stabilité et FPS.
- Sur mobile/ultrabook, réduire **Qualité** pour fluidité.

---

## 5) Conseils d’utilisation (quick start)

1. Commencer en mode `NORMAL`.
2. Régler `Inclinaison` et `Rayon ext.`.
3. Ajuster `Lensing` et `Spin a` progressivement.
4. Activer `BLOOM` pour un rendu plus cinématique.
5. Si FPS bas : baisser `Qualité` et/ou `Rayon ext.`.

---

## 6) Limites actuelles (connues)

- Modèle physique volontairement "Kerr-lite" / stylisé (objectif visuel temps réel).
- Pas de pipeline build/test automatisé (projet statique simple).
- Pas de persistance des presets utilisateur dans cette version.

---

## 7) Roadmap suggérée (sans casser l’existant)

- Ajouter système de presets (JSON localStorage).
- Ajouter capture d’écran intégrée.
- Ajouter panneau "aide" dépliable.
- Ajouter un mode colorimétrique scientifique vs cinématique.

---

## 8) Fichiers du projet

- `index.html`
- `style.css`
- `script.js`
- `README.md`

---

## 9) Licence

Aucune licence explicite détectée dans le dépôt.

Si ce projet doit être partagé publiquement, ajouter un fichier `LICENSE` (MIT, Apache-2.0, etc.) selon ton intention.
