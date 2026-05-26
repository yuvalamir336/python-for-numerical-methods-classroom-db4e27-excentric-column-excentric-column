def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """
    import numpy as np
from scipy.optimize import bisect
   
    # 1. חישוב עומס הקריסה התיאורטי של אוילר כחסם עליון
    P_euler = (np.pi**2 * E * A) / (L / r) ** 2

    # 2. הגדרת פונקציית המטרה למציאת שורש: f(P) = sigma_max(P) - sigma_allow
    def f(P):
        if P == 0:
            return -sigma_allow

        # הארגומנט של הקוסינוס (ברדיאנים)
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # חישוב ערך הסקנט (ההופכי של קוסינוס)
        sec_val = 1.0 / np.cos(theta)

        # נוסחת הסקנט למאמץ המקסימלי בסיב הקיצוני
        sigma_max = (P / A) * (1 + (e * c / r**2) * sec_val)

        return sigma_max - sigma_allow

    # 3. פתרון נומרי בשיטת החצייה (Bisection)
    # החסם העליון נקבע מעט מתחת ל-P_euler כדי למנוע חלוקה באפס בתוך הקוסינוס
    P_critical = bisect(f, 0, 0.9999 * P_euler)

    return float(P_critical)


# ==============================================================================
# חלק הרצה - כאן את יכולה להזין את הנתונים שלך כדי שהתוכנה תריץ את החישוב:
# ==============================================================================
if __name__ == "__main__":
    # דוגמה לערכים (שני את המספרים האלה בהתאם לפרופיל שלך):
    L_input = 3000.0  # אורך המוט במ"מ
    E_input = 200000.0  # מודול אלסטיות ב-MPa (עבור פלדה)
    A_input = 5000.0  # שטח חתך בממ"ר
    r_input = 50.0  # רדיוס אינרציה במ"מ
    c_input = 100.0  # מרחק לסיב קיצוני במ"מ
    e_input = 20.0  # אקסצנטריות העומס במ"מ
    sigma_allow_input = 250.0  # מאמץ מותר ב-MPa

    # הרצת הפונקציה
    try:
        result_P = find_critical_load(
            L=L_input,
            E=E_input,
            A=A_input,
            r=r_input,
            c=c_input,
            e=e_input,
            sigma_allow=sigma_allow_input,
        )

        print("-" * 50)
        print(f"The Critical Axial Load (P) is: {result_P:.2f} Newtons")
        print(f"In Kilo-Newtons (kN): {result_P / 1000.0:.2f} kN")
        print("-" * 50)

    except ValueError as err:
        print(f"Numerical convergence error: {err}")
