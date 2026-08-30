from __future__ import annotations

import pytest
from PySide6.QtWidgets import QApplication, QDialog

from ddmsoft.fitting import get_model
from ddmsoft.gui.fit_dialog import InitialGuessDialog


@pytest.fixture
def qapp():
    instance = QApplication.instance() or QApplication([])
    yield instance
    for widget in instance.topLevelWidgets():
        widget.close()


def test_initial_guess_dialog_is_generated_from_model_and_validates_values(qapp):
    model = get_model("stretch")
    dialog = InitialGuessDialog(model, model.defaults, (False,) * len(model.parameters))

    assert len(dialog._value_edits) == len(model.parameters)
    assert dialog._value_edits[0].text() == "2e-12"
    assert dialog._value_edits[-1].text() == ""

    dialog._value_edits[0].setText("not numeric")
    dialog._validate_and_accept()
    assert dialog.result() == QDialog.DialogCode.Rejected
    assert "must be numeric" in dialog.error_label.text()

    dialog._value_edits[0].setText("3e-12")
    dialog._value_edits[-1].setText("0.02")
    dialog._fixed_checks[-1].setChecked(True)
    dialog._validate_and_accept()

    assert dialog.result() == QDialog.DialogCode.Accepted
    assert dialog.initial_values == (3e-12, 1.0, None, 0.02)
    assert dialog.fixed_flags == (False, False, False, True)


def test_initial_guess_dialog_rejects_fixed_automatic_values(qapp):
    model = get_model("cumulant_1")
    dialog = InitialGuessDialog(model, model.defaults, (False,) * len(model.parameters))
    dialog._fixed_checks[-1].setChecked(True)
    dialog._validate_and_accept()

    assert dialog.result() == QDialog.DialogCode.Rejected
    assert "needs a numeric value" in dialog.error_label.text()
