# JAR LINE Coding Style Guide

## Code Style and Conventions

### Infinity and Maximum Value Checks
- **Replace `== Double.POSITIVE_INFINITY` with `Utils.isInf()`**
- **Replace `== Integer.MAX_VALUE` with `Utils.isInf()`**
- Always use the utility method instead of direct comparisons to maintain consistency

### Documentation Requirements
- **Always add Javadoc comments** for all public methods, classes, and fields
- Javadoc generated with assistance from a LLM or similar tools is acceptable

### Protected Files and Components
- **Do NOT add new methods to `Matrix.java` without explicit approval**
  - Always ask before modifying this core class
  - Discuss the necessity and design of any proposed additions
- **Do NOT modify `pom.xml` without permission**
  - Maven configuration changes require review
  - Ask before adding dependencies or changing build settings

### Code Quality Standards
- **Do not add unnecessary code or temporary tests to the codebase**

### Migration and Incomplete Work
- **Always leave TODO comments when migration is incomplete**
- Format: `// TODO: [Description of what needs to be completed]`

### Testing Requirements
- **Add unit tests as you develop** - don't leave testing for later
- **Use JUnit 5 exclusively**
- **Import statement: `import org.junit.jupiter.api.Test;`**
- Do not use other testing frameworks or older JUnit versions
- Write clear, descriptive test method names

### Logging and Output Standards
- **Do NOT use `System.out.print*` or `System.err.print*` methods**
- **Always use LINE's logging functions instead:**
  - `line_printf(format, args...)` for debug/info output
  - `line_warning(mfilename(new Object() {}), message)` for warnings
  - `line_error(mfilename(new Object() {}), message)` for errors
- **Use `mfilename(new Object() {})` for proper method identification in warnings/errors**
- This ensures consistent formatting and integration with LINE's verbose level settings

