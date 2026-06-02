import os
from PyQt5.QtWidgets import (
    QWidget, QGridLayout, QVBoxLayout, QLabel,
    QComboBox, QGroupBox, QSizePolicy,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap

_HERE = os.path.dirname(os.path.abspath(__file__))
PNG_PATH = os.path.join(_HERE, 'circuit.png')

# ── Capacitor value tables ─────────────────────────────────────────────────────
ALL_DISPLAY = ['1 pF', '10 pF', '100 pF', '1000 pF', '10 nF', '100 nF', '1 µF', '10 µF']
ALL_CFGKEYS = ['1pF',  '10pF',  '100pF',  '1000pF',  '10nF',  '100nF',  '1uF',  '10uF']

TOP_DISPLAY = ALL_DISPLAY[1:]   # 7 items: 10pF … 10µF
TOP_CFGKEYS = ALL_CFGKEYS[1:]

BOT_DISPLAY = ALL_DISPLAY[:-1]  # 7 items: 1pF … 1µF
BOT_CFGKEYS = ALL_CFGKEYS[:-1]

# Constraint: TOP_CFGKEYS[k] = 10 × BOT_CFGKEYS[k], so both combos share the same index.

# ── Capacitor categories ───────────────────────────────────────────────────────
CATEGORIES = ['AIR', 'CERAMIC', 'HIGHVALUE']

# Physical value set for each category (includes sub-measurable values for classification)
CATEGORY_VALUES = {
    'AIR':       {'1pF', '10pF', '100pF', '1000pF'},
    'CERAMIC':   {'10nF', '100nF', '1uF'},
    'HIGHVALUE': {'10uF'},
}
VALUE_TO_CATEGORY = {v: cat for cat, vals in CATEGORY_VALUES.items() for v in vals}


def _valid_indices(cfgkeys: list, category) -> set:
    """Indices in cfgkeys whose value belongs to category (None = all enabled)."""
    if category is None:
        return set(range(len(cfgkeys)))
    allowed = CATEGORY_VALUES.get(category, set())
    return {i for i, k in enumerate(cfgkeys) if k in allowed}


# ── PNG canvas ─────────────────────────────────────────────────────────────────
class PDFCanvas(QLabel):
    """Displays a PNG image, preserving aspect ratio on resize."""

    def __init__(self, png_path: str):
        super().__init__()
        self.setAlignment(Qt.AlignCenter)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setMinimumSize(200, 200)
        self._source = QPixmap(png_path)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not self._source.isNull():
            self.setPixmap(
                self._source.scaled(self.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )


# ── Setup widget ───────────────────────────────────────────────────────────────
class CircuitSetupWidget(QWidget):
    """Setup tab: bridge schematic with capacitance and serial-number selectors."""

    def __init__(self, config_path: str):
        super().__init__()
        self.config_path = config_path
        self._sn_categories = {}                    # SN (upper) -> category name
        self._cat_serials = {c: [] for c in CATEGORIES}  # category -> [SNs]
        self._build_ui()
        self._load()

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)

        mid = QWidget()
        grid = QGridLayout(mid)
        grid.setSpacing(10)

        self.schematic = PDFCanvas(PNG_PATH)
        grid.addWidget(self.schematic, 0, 0, 2, 1)

        self._gb1, self.cb1_val, self.cb1_sn = self._cap_group('C1', top=True)
        self._gb2, self.cb2_val, self.cb2_sn = self._cap_group('C2', top=False)
        grid.addWidget(self._gb1, 0, 1, Qt.AlignTop)
        grid.addWidget(self._gb2, 1, 1, Qt.AlignBottom)

        outer.addWidget(mid)

        # 10:1 constraint — skip sync if target index is disabled by category filter
        self.cb1_val.currentIndexChanged.connect(lambda i: self._sync(self.cb2_val, i))
        self.cb2_val.currentIndexChanged.connect(lambda i: self._sync(self.cb1_val, i))

        # Category filter updates when SN changes
        self.cb1_sn.currentTextChanged.connect(lambda t: self._update_filter(self.cb1_val, t, True))
        self.cb2_sn.currentTextChanged.connect(lambda t: self._update_filter(self.cb2_val, t, False))

    def _cap_group(self, title: str, top: bool):
        box = QGroupBox(title)
        lay = QVBoxLayout(box)
        val_combo = QComboBox()
        for lbl in (TOP_DISPLAY if top else BOT_DISPLAY):
            val_combo.addItem(lbl)
        sn_combo = QComboBox()
        sn_combo.setEditable(True)
        sn_combo.setInsertPolicy(QComboBox.NoInsert)
        sn_combo.setMinimumWidth(130)
        lay.addWidget(QLabel('Capacitance:'))
        lay.addWidget(val_combo)
        lay.addWidget(QLabel('Serial no.:'))
        lay.addWidget(sn_combo)
        return box, val_combo, sn_combo

    def _sync(self, other: QComboBox, idx: int):
        """Sync 10:1 partner — only if that index is enabled for the other SN's category."""
        item = other.model().item(idx)
        if item and (item.flags() & Qt.ItemIsEnabled):
            other.blockSignals(True)
            other.setCurrentIndex(idx)
            other.blockSignals(False)

    def _update_filter(self, val_combo: QComboBox, sn_text: str, is_top: bool):
        """Enable only values belonging to the selected SN's category."""
        category = self._sn_categories.get(sn_text.strip().upper())
        cfgkeys = TOP_CFGKEYS if is_top else BOT_CFGKEYS
        valid = _valid_indices(cfgkeys, category)

        model = val_combo.model()
        for i in range(val_combo.count()):
            item = model.item(i)
            if i in valid:
                item.setFlags(item.flags() | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            else:
                item.setFlags(item.flags() & ~(Qt.ItemIsEnabled | Qt.ItemIsSelectable))

        # If the current selection is now disabled, jump to the first valid index
        if val_combo.currentIndex() not in valid and valid:
            val_combo.blockSignals(True)
            val_combo.setCurrentIndex(min(valid))
            val_combo.blockSignals(False)

    # ── Config I/O ─────────────────────────────────────────────────────────────

    def _parse_raw_config(self):
        """Parse config.ini.

        Returns (all_sns, active, sn_categories, cat_serials).
        all_sns  — ordered list of every known SN (category sections first, then [CONFIG])
        active   — dict of active key→value from [CONFIG]
        sn_categories — SN (upper) → category name
        cat_serials   — category → [SNs] from category sections
        """
        active = {}
        sn_categories = {}
        cat_serials = {c: [] for c in CATEGORIES}
        config_sns = []     # SNs seen in [CONFIG] (active + commented)
        current_section = None

        with open(self.config_path, 'r') as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith('[') and stripped.endswith(']'):
                    current_section = stripped[1:-1].upper()
                    continue

                if current_section == 'CONFIG':
                    if not stripped:
                        continue
                    is_comment = stripped.startswith(('#', ';'))
                    content = stripped.lstrip('#; ').strip()
                    if '=' not in content:
                        continue
                    key, val = content.split('=', 1)
                    key, val = key.strip().lower(), val.strip()
                    if not val:
                        continue
                    if key.startswith('sn'):
                        val_up = val.upper()
                        if val_up not in config_sns:
                            config_sns.append(val_up)
                        if not is_comment:
                            active[key] = val_up
                    elif not is_comment:
                        active[key] = val

                elif current_section in CATEGORIES:
                    if stripped and not stripped.startswith(('#', ';')) and '=' in stripped:
                        key, val = stripped.split('=', 1)
                        if key.strip().lower() == 'serials':
                            for sn in val.split(','):
                                sn = sn.strip().upper()
                                if sn:
                                    if sn not in cat_serials[current_section]:
                                        cat_serials[current_section].append(sn)
                                    sn_categories[sn] = current_section

        # Ordered SN list: category sections first, then any from [CONFIG] not yet seen
        seen = set()
        all_sns = []
        for cat in CATEGORIES:
            for sn in cat_serials[cat]:
                if sn not in seen:
                    all_sns.append(sn)
                    seen.add(sn)
        for sn in config_sns:
            if sn not in seen:
                all_sns.append(sn)
                seen.add(sn)

        return all_sns, active, sn_categories, cat_serials

    def _load(self):
        all_sns, active, sn_categories, cat_serials = self._parse_raw_config()
        self._sn_categories = sn_categories
        self._cat_serials = cat_serials

        for cb in (self.cb1_sn, self.cb2_sn):
            cb.blockSignals(True)
            cb.clear()
            for sn in all_sns:
                cb.addItem(sn)
            cb.blockSignals(False)

        # Capacitance values (set before SNs so filter has the right starting point)
        for cfg_key, combo, is_top in [
            ('c1',    self.cb1_val, True),
            ('c2', self.cb2_val, False),
        ]:
            raw = active.get(cfg_key, '').upper().replace(' ', '')
            keys = [k.upper() for k in (TOP_CFGKEYS if is_top else BOT_CFGKEYS)]
            try:
                combo.blockSignals(True)
                combo.setCurrentIndex(keys.index(raw))
                combo.blockSignals(False)
            except ValueError:
                pass

        # Serial numbers + apply initial category filter
        for cfg_key, sn_combo, val_combo, is_top in [
            ('snc1',   self.cb1_sn, self.cb1_val, True),
            ('snc2',self.cb2_sn, self.cb2_val, False),
        ]:
            sn = active.get(cfg_key, '')
            idx = sn_combo.findText(sn, Qt.MatchFixedString)
            if idx >= 0:
                sn_combo.blockSignals(True)
                sn_combo.setCurrentIndex(idx)
                sn_combo.blockSignals(False)
            self._update_filter(val_combo, sn, is_top)

    def save_config(self):
        """Write selections to [CONFIG] and assign any new SNs to their category section."""
        positions = [
            ('c1',    'snc1',   self.cb1_val, self.cb1_sn, True),
            ('c2', 'snc2',self.cb2_val, self.cb2_sn, False),
        ]

        updates = {}
        new_assignments = {}   # SN -> category (only for truly new SNs)

        for val_key, sn_key, val_combo, sn_combo, is_top in positions:
            cfgkeys = TOP_CFGKEYS if is_top else BOT_CFGKEYS
            val = cfgkeys[val_combo.currentIndex()]
            sn = sn_combo.currentText().strip().upper()
            updates[val_key] = val
            updates[sn_key] = sn
            if sn and sn not in self._sn_categories:
                cat = VALUE_TO_CATEGORY.get(val)
                if cat:
                    new_assignments[sn] = cat

        # Update [CONFIG] section, preserving comments
        with open(self.config_path, 'r') as f:
            lines = f.readlines()

        in_config = False
        config_end = len(lines)
        updated = set()

        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.lower() == '[config]':
                in_config = True
                continue
            if stripped.startswith('[') and in_config:
                config_end = i
                in_config = False
                break
            if in_config and not stripped.startswith(('#', ';')) and '=' in stripped:
                key = stripped.split('=', 1)[0].strip().lower()
                if key in updates:
                    lines[i] = f'{key} = {updates[key]}\n'
                    updated.add(key)

        for key, val in updates.items():
            if key not in updated:
                lines.insert(config_end, f'{key} = {val}\n')
                config_end += 1

        with open(self.config_path, 'w') as f:
            f.writelines(lines)

        # Assign new SNs to their category section
        if new_assignments:
            for sn, cat in new_assignments.items():
                if sn not in self._cat_serials[cat]:
                    self._cat_serials[cat].append(sn)
                self._sn_categories[sn] = cat
            self._update_category_sections()

    def _update_category_sections(self):
        """Rewrite the serials line in each [AIR]/[CERAMIC]/[HIGHVALUE] section."""
        with open(self.config_path, 'r') as f:
            lines = f.readlines()

        for cat in CATEGORIES:
            header_lower = f'[{cat.lower()}]'
            serials_line = 'serials = ' + ', '.join(self._cat_serials[cat]) + '\n'

            section_idx = next(
                (i for i, l in enumerate(lines) if l.strip().lower() == header_lower),
                None
            )

            if section_idx is None:
                lines.append(f'\n[{cat}]\n')
                lines.append(serials_line)
            else:
                i = section_idx + 1
                found = False
                while i < len(lines):
                    s = lines[i].strip()
                    if s.startswith('['):
                        break
                    if s and not s.startswith(('#', ';')) and '=' in s:
                        if s.split('=', 1)[0].strip().lower() == 'serials':
                            lines[i] = serials_line
                            found = True
                            break
                    i += 1
                if not found:
                    lines.insert(i, serials_line)

        with open(self.config_path, 'w') as f:
            f.writelines(lines)
