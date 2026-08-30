"""Dialogs for collecting registry-defined fit parameters."""

from __future__ import annotations

from collections.abc import Sequence
from math import isfinite

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
)

from ..fitting import ModelDefinition


class InitialGuessDialog(QDialog):
    """Edit initial values and fixed flags for one registered model."""

    def __init__(
        self,
        model: ModelDefinition,
        initial_values: Sequence[float | None],
        fixed_flags: Sequence[bool],
        *,
        parent=None,
    ) -> None:
        super().__init__(parent)
        values = tuple(initial_values)
        fixed = tuple(bool(value) for value in fixed_flags)
        if len(values) != len(model.parameters) or len(fixed) != len(model.parameters):
            raise ValueError("initial guess state does not match the selected model")
        self.model = model
        self._value_edits: list[QLineEdit] = []
        self._fixed_checks: list[QCheckBox] = []
        self.initial_values: tuple[float | None, ...] = values
        self.fixed_flags: tuple[bool, ...] = fixed
        self.setWindowTitle(f"Initial guess: {model.display_name}")
        self.setModal(True)

        layout = QVBoxLayout(self)
        description = QLabel(
            "Enter finite numeric values. Leave amplitude or background blank to estimate it."
        )
        description.setWordWrap(True)
        layout.addWidget(description)
        form = QFormLayout()
        for index, (definition, value, is_fixed) in enumerate(zip(model.parameters, values, fixed)):
            edit = QLineEdit()
            edit.setObjectName(f"initialValueEdit_{definition.identifier}")
            edit.setPlaceholderText("automatic estimate" if value is None else "value")
            if value is not None:
                edit.setText(f"{value:.12g}")
            check = QCheckBox("Fixed")
            check.setObjectName(f"fixedCheck_{definition.identifier}")
            check.setChecked(is_fixed)
            row = QHBoxLayout()
            row.addWidget(edit, 1)
            row.addWidget(check)
            form.addRow(definition.display_name, row)
            self._value_edits.append(edit)
            self._fixed_checks.append(check)
        layout.addLayout(form)

        self.error_label = QLabel()
        self.error_label.setObjectName("initialGuessError")
        self.error_label.setWordWrap(True)
        self.error_label.setStyleSheet("color: #b00020")
        layout.addWidget(self.error_label)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _validate_and_accept(self) -> None:
        values: list[float | None] = []
        for definition, edit, check in zip(
            self.model.parameters, self._value_edits, self._fixed_checks
        ):
            text = edit.text().strip()
            if not text:
                if check.isChecked():
                    self.error_label.setText(
                        f"{definition.display_name}: a fixed parameter needs a numeric value."
                    )
                    return
                values.append(None)
                continue
            try:
                value = float(text)
            except ValueError:
                self.error_label.setText(f"{definition.display_name}: value must be numeric.")
                return
            if not isfinite(value):
                self.error_label.setText(f"{definition.display_name}: value must be finite.")
                return
            values.append(value)
        self.initial_values = tuple(values)
        self.fixed_flags = tuple(check.isChecked() for check in self._fixed_checks)
        self.error_label.clear()
        self.accept()


__all__ = ["InitialGuessDialog"]
