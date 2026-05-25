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
    # חישוב עומס הקריסה של אוילר
    # Pcr = (π²EA)/(L/r)²
    P_euler = (np.pi*2 * E * A) / (L / r)*2

    # פונקציית העזר:
    # f(P) = σmax(P) - σallow
    def f(P):

        # במקרה של עומס אפס
        if P == 0:
            return -sigma_allow

        # הארגומנט של הקוסינוס (ברדיאנים)
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))

        # sec(theta) = 1 / cos(theta)
        sec_val = 1.0 / np.cos(theta)

        # נוסחת הסקנט למאמץ מקסימלי
        sigma_max = (
            (P / A)
            * (1 + (e * c / r**2) * sec_val)
        )

        # מחפשים f(P)=0
        return sigma_max - sigma_allow

    # פתרון נומרי בשיטת החצייה
    # נמנעים מהאסימפטוטה ב-P=P_euler
    P_critical = bisect(
        f,
        0,
        0.9999 * P_euler
    )

    return float(P_critical)
